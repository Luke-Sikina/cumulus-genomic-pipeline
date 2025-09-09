import argparse
import logging
from pathlib import Path
from pydantic import BaseModel

class VcfProcessingInput(BaseModel):
    vcf_files: list[str]
    output_dir: str
    valid: bool

def validate(args: argparse.Namespace) -> VcfProcessingInput:
    logging.info("Validating CLI args...")
    vcf_files: list[str] = []
    if 'vcf' in args and args.vcf:
        logging.info(f'Found {len(args.vcf)} vcf args')
        for vcf_file in args.vcf:
            vcf_path = Path(vcf_file)
            if not vcf_path.exists():
                logging.warning(f"{vcf_file} does not exist")
            elif not vcf_path.is_file():
                logging.warning(f"{vcf_file} is not a regular file")
            else:
                logging.info(f"VCF {vcf_file} exists and is normal file")
                vcf_files.append(vcf_file)
        logging.info(f"{len(vcf_files)} valid input files found.")
    else:
        logging.error('No vcf arg found')
    valid = len(vcf_files) > 0

    output_dir = args.output_dir if 'output_dir' in args else ''
    logging.info(f'Validating output dir {output_dir}')
    if output_dir:
        logging.info(f"Output configured to {output_dir}")
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        if not output_path.exists():
            logging.error(f"Output dir {output_dir} does not exist and cannot be created.")
            valid = False
        logging.info(f"Confirmed output dir {output_dir} exists or was created")
    else:
        logging.error("No output_dir argument found")
        valid = False
    
    return VcfProcessingInput(vcf_files=vcf_files, output_dir=output_dir, valid=valid)

