#!/usr/bin/env python3
import argparse
import sys
import logging
from pathlib import Path

from cumulus_genomic_pipeline.process_args import validate
from cumulus_genomic_pipeline.process_vcf import process_inputs

def main():
    parser = argparse.ArgumentParser(description="Process multiple files with Docker")    
    parser.add_argument('-i', '--vcf', action='append', required=True,
                       help='Input file paths (specify multiple times)')
    parser.add_argument('-o', '--output_dir', required=True,
                       help='Output directory path')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output')
    args = parser.parse_args()
    
    if args.verbose:
        logging.basicConfig(level=logging.INFO)
    
    inputs = validate(args)
    if inputs.valid:
        logging.info('Valid CLI args. Processing VCFs...')
        process_inputs(inputs)
        logging.info('Done processing. Exiting...')
        sys.exit(0)
    else:
        logging.error('Invalid CLI args. Exiting...')
        sys.exit(1)

if __name__ == "__main__":
    main()