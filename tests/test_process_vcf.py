import argparse
from pathlib import PosixPath

from cumulus_genomic_pipeline.schema.schema import variant_schema, consequence_schema, occurance_schema
from cumulus_genomic_pipeline.process_args import VcfProcessingInput, validate
from cumulus_genomic_pipeline.process_vcf import CONSEQUENCE_OUT, OCCURANCE_OUT, VARIANT_OUT, process_inputs
from tests.utils.utils import verify_parquet_file


def test_process_vcf(tmp_path):
    # Simple e2e test. Does the program read in a VCF, process it, and
    # output to a parquet file?
    output_dir: PosixPath = tmp_path / "output"
    output_dir.mkdir()

    inputs = VcfProcessingInput(
        vcf_files=['tests/data/4klines.variants.CEPH-1463.snv.vep.vcf.gz'],
        output_dir=f"{output_dir.resolve()}",
        valid=False
    )
    
    process_inputs(inputs)
    
    variant = output_dir / VARIANT_OUT
    consequence = output_dir / CONSEQUENCE_OUT
    occurance = output_dir / OCCURANCE_OUT
    assert variant.exists()
    assert consequence.exists()
    assert occurance.exists()
    
    # These counts are more assertions that there are a non-zero number of rows
    # than they are assertions that there should be that exact no. of rows
    assert verify_parquet_file(variant, variant_schema, 561)
    assert verify_parquet_file(consequence, consequence_schema, 4443)
    assert verify_parquet_file(occurance, occurance_schema, 561)
