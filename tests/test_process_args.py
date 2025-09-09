import argparse
from pathlib import PosixPath

from cumulus_genomic_pipeline.process_args import VcfProcessingInput, validate


def test_fail_empty_args():
    args = argparse.Namespace()
    actual: VcfProcessingInput = validate(args)
    expected = VcfProcessingInput(vcf_files=[], output_dir="", valid=False)

    assert expected == actual

def test_fail_no_vcf_files(tmp_path):
    output_dir: PosixPath = tmp_path / "output"
    output_dir.mkdir()
    
    args = argparse.Namespace(
        vcf=[],
        output_dir=f"{output_dir.resolve()}"
    )
    
    actual: VcfProcessingInput = validate(args)
    expected = VcfProcessingInput(vcf_files=[], output_dir=f"{output_dir.resolve()}", valid=False)
    
    assert actual == expected

def test_fail_no_out_dir(tmp_path):
    input_dir: PosixPath = tmp_path / "input"
    input_dir.mkdir()
    vcf_file_a = input_dir / "sample1.vcf"
    vcf_file_b = input_dir / "sample2.vcf"
        
    vcf_file_a.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n")    
    args = argparse.Namespace(
        vcf=[f"{vcf_file_a}", f"{vcf_file_b}"],
    )
    
    actual: VcfProcessingInput = validate(args)
    #                                          2nd vcf filtered b/c DNE
    expected = VcfProcessingInput(vcf_files=[f"{vcf_file_a}"], output_dir='', valid=False)
    
    assert actual == expected
