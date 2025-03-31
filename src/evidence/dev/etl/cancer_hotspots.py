"""Module for ETL cancer hotspots data"""

import datetime
import json
import logging
from pathlib import Path
from timeit import default_timer as timer
from typing import AsyncGenerator, Dict

import pandas as pd
import requests
# from variation.query import QueryHandler

from evidence import DATA_DIR_PATH
from evidence.data_sources import CancerHotspots

from ga4gh.core.models import MappableConcept, Coding
from ga4gh.vrs.models import Allele, Location
from ga4gh.cat_vrs.recipes import ProteinSequenceConsequence
from ga4gh.cat_vrs.models import CategoricalVariant, DefiningLocationConstraint, DefiningAlleleConstraint
from ga4gh.va_spec.base import CohortAlleleFrequencyStudyResult

class MockQueryHandlerResponse:
    def __init__(self, variation):
        self.variation=variation

class MockQueryHandler:
    async def normalize(self, variation: str):
        from ga4gh.vrs.models import SequenceLocation, LiteralSequenceExpression
        
        variation = Allele(
            location=SequenceLocation(
                start=1,
                end=1,
                sequenceReference="NP_001191.1"
            ),
            state=LiteralSequenceExpression(sequence="A")
        )
        return MockQueryHandlerResponse(variation=variation)

    def __init__(self):
        self.normalize_handler = self

class CancerHotspotsETLError(Exception):
    """Exceptions for Cancer Hotspots ETL"""


_logger = logging.getLogger(__name__)


class CancerHotspotsETL(CancerHotspots):
    """Class for Cancer Hotspots ETL methods."""

    _cat_var_relations = [
        MappableConcept(
            primaryCoding=Coding(
                code="translation_of",
                system="http://www.sequenceontology.org"
            )
        )
    ]
    _cat_var_match_characteristic = MappableConcept(
        primaryCoding=Coding(
            code="exactMatch",
            system="http://www.sequenceontology.org" # TODO
        )
    )

    def __init__(
        self,
        data_url: str = "https://www.cancerhotspots.org/files/hotspots_v2.xls",
        src_dir_path: Path = DATA_DIR_PATH / "cancer_hotspots",
        transformed_data_path: Path | None = DATA_DIR_PATH / "cancer_hotspots",
        ignore_transformed_data: bool = True,
    ) -> None:
        """Initialize the CancerHotspotsETL class

        :param data_url: URL to data file
        :param src_dir_path: Path to cancer hotspots data directory
        :param transformed_data_path: Path to transformed cancer hotspots file
        :param ignore_transformed_data: `True` if only bare init is needed. This is
            intended for developers when using the CLI to transform cancer hotspots
            data. Ignores path set in `_transformed_data_path`. `False` will load
            transformed data from s3
        """
        super().__init__(
            data_url, src_dir_path, transformed_data_path, ignore_transformed_data
        )
        fn = self.data_url.split("/")[-1]
        self.data_path = self.src_dir_path / fn

    def download_data(self) -> None:
        """Download Cancer Hotspots data."""
        if not self.data_path.exists():
            r = requests.get(self.data_url, timeout=5)
            if r.status_code == 200:
                with self.data_path.open("wb") as f:
                    f.write(r.content)
            else:
                _logger.error(
                    "Unable to download Cancer Hotspots data. Received status code: %i",
                    r.status_code,
                )

    async def transform_and_write_hotspots(self) -> None:
        """Normalize variations in cancer hotspots and updates `transformed_data`

        Run manually each time variation-normalizer or Cancer Hotspots releases a new
        version.
        """
        self.download_data()
        if not self.data_path.exists():
            err_msg = "Downloading Cancer Hotspots data was unsuccessful"
            raise CancerHotspotsETLError(err_msg)

        snv_hotspots = pd.read_excel(self.data_path, sheet_name="SNV-hotspots")
        indel_hotspots = pd.read_excel(self.data_path, sheet_name="INDEL-hotspots")
        variation_normalizer = MockQueryHandler() # QueryHandler()

        _logger.info("Normalizing Cancer Hotspots data...")

        today = datetime.datetime.strftime(
            datetime.datetime.now(tz=datetime.UTC), "%Y%m%d"
        )
        transformed_data_path = self.src_dir_path / f"cancer_hotspots_{today}.json"
        with transformed_data_path.open("w") as f:
            start = timer()
            async for study_result in self.get_transformed_data(snv_hotspots, variation_normalizer, is_snv=True):
                f.write(study_result.model_dump_json())
                f.write("\n")

            async for study_result in self.get_transformed_data(indel_hotspots, variation_normalizer, is_snv=False):
                f.write(study_result.model_dump_json())
                f.write("\n")

        end = timer()
        _logger.info("Successfully transformed Cancer Hotspots data in %.*f s", 2, end - start)


    async def get_transformed_data(
        self, df: pd.DataFrame, variation_normalizer, is_snv: bool
    ) -> AsyncGenerator[CohortAlleleFrequencyStudyResult, None]:
        """Normalize variant and updates `transformed_data`

        :param df: Dataframe to transform
        :param variation_normalizer: Variation Normalizer handler
        :param is_snv: `True` if SNV data, else INDEL
        """
        aa_location_group_cols = ["Hugo_Symbol", "Amino_Acid_Position"]
        grouped_df = df.groupby(aa_location_group_cols).apply(lambda x: x.to_dict('records'), include_groups=False).reset_index(name='Rows')
        for _, group in grouped_df.iterrows(): 
            normalized_allele = None
            for row in group["Rows"]:
                normalized_allele = await self.normalize_row(variation_normalizer, is_snv, group, row)
                prot_cons_cat_var = self._create_protein_seq_cons(normalized_allele)
                row["ProteinSequenceConsequence"] = prot_cons_cat_var

            def_loc_cat_var = self._create_loc_cat_var(normalized_allele.location)
            group["DefiningLocationCatVar"] = def_loc_cat_var

            yield self._create_study_result(group)
            
    def _create_loc_cat_var(self, loc: Location):
        def_loc_constraint = DefiningLocationConstraint(
            location=loc,
            relations=CancerHotspotsETL._cat_var_relations,
            matchCharacteristic=CancerHotspotsETL._cat_var_match_characteristic
        )
        return CategoricalVariant(constraints=[def_loc_constraint])
    
    def _create_protein_seq_cons(self, p_allele: Allele):
        def_allele_constraint = DefiningAlleleConstraint(
            allele=p_allele,
            relations=CancerHotspotsETL._cat_var_relations
        )
        return ProteinSequenceConsequence(constraints=[def_allele_constraint])

    def _create_study_result(self, group: Dict):
        return group["DefiningLocationCatVar"]

    async def normalize_row(self, variation_normalizer, is_snv: bool, group: Dict, row: Dict):
        hugo_symbol = group["Hugo_Symbol"]
        pos = group["Amino_Acid_Position"]
        alt = row["Variant_Amino_Acid"].split(':')[0]

        if is_snv:
            ref = row["ref"]
            variation = f"{hugo_symbol} {ref}{pos}{alt.split(':')[0]}"
        else:
            ref = None
            variation = f"{hugo_symbol} {alt.split(':')[0]}"

        try:
            variation_norm_resp = (
                    await variation_normalizer.normalize_handler.normalize(variation)
                )

        except Exception as e:
            _logger.error(
                    "variation-normalizer unable to normalize %s: %s", variation, str(e)
                )

        else:
            if variation_norm_resp and variation_norm_resp.variation:
                return variation_norm_resp.variation
                
            else:
                _logger.warning(
                        "variation-normalizer unable to normalize: %s", variation
                    )
    
