import logging

import pyarrow as pa
import pyarrow.parquet as pq
from pyarrow import Schema
from pathlib import Path

from cyvcf2 import VCF, Variant
from cumulus_genomic_pipeline.radiant.vcf.experiment import Case, Experiment
from cumulus_genomic_pipeline.schema.schema import variant_schema, consequence_schema, occurance_schema
from cumulus_genomic_pipeline.process_args import VcfProcessingInput
from cumulus_genomic_pipeline.radiant.vcf.common import process_common
from cumulus_genomic_pipeline.radiant.vcf.consequence import parse_csq_header, process_consequence
from cumulus_genomic_pipeline.radiant.vcf.occurrence import process_occurrence
from cumulus_genomic_pipeline.radiant.vcf.pedigree import Pedigree
from cumulus_genomic_pipeline.radiant.vcf.variant import process_variant

VARIANT_OUT = 'variants.parquet'
OCCURANCE_OUT = 'occurance.parquet'
CONSEQUENCE_OUT = 'consequence.parquet'
BATCH_SIZE = 1000

def process_inputs(inputs: VcfProcessingInput):
    case_id = 0
    for vcf_path in inputs.vcf_files:
        case_id+=1
        _process_vcf(vcf_path, inputs.output_dir, case_id)

def _process_vcf(vcf_path: str, output_dir: str, case_id: int):
    variant_pq = Path(output_dir) / VARIANT_OUT
    occurance_pq = Path(output_dir) / OCCURANCE_OUT
    consequence_pq = Path(output_dir) / CONSEQUENCE_OUT
    logging.info(f"Processing vcf {vcf_path} outputting to {output_dir}")
    
    with pq.ParquetWriter(variant_pq, variant_schema) as variant_writer, \
        pq.ParquetWriter(consequence_pq, consequence_schema) as conseq_writer, \
        pq.ParquetWriter(occurance_pq, occurance_schema) as occurance_writer:
        vcf = VCF(vcf_path)
        logging.debug(f"Cases: {vcf.samples}")
        csq_header = parse_csq_header(vcf)
        experiments: list[Experiment] = []
        for sample in vcf.samples:
            experiments.append(Experiment(seq_id=1, task_id=1, patient_id=1, aliquot=sample, family_role='child', affected_status='', sex='Unknown', experimental_strategy='Unknown'))
        case = Case(case_id=1, part=1, vcf_filepath=vcf_path, analysis_type='WGS', experiments=experiments, index_vcf_filepath=None)
        ped = Pedigree(case, vcf.samples)
        logging.info(f'Found the following samples: {vcf.samples}')
        record_count = 0
        record: Variant
        variant_batch = []
        occurance_batch = []
        consequence_batch = []
        for record in vcf:
            record_count += 1
            if record_count % 1000 == 0:
                logging.debug(f'Record count: {record_count} vars: {len(variant_batch)} cons: {len(consequence_batch)} occ: {len(occurance_batch)}')
                _write_all_tables(variant_writer, conseq_writer, occurance_writer, variant_batch, consequence_batch, occurance_batch)
                variant_batch = []
                occurance_batch = []
                consequence_batch = []
            out =_process_record(case_id, csq_header, ped, record, vcf_path)
            if out is None:
                logging.warning('Discarding record #{record_count}')
            else:
                (consequences, occurrences, variant) = out
                variant_batch.append(variant)
                consequence_batch.extend(consequences)
                occurance_batch.extend(occurrences.values())
        logging.debug(f'Record count: {record_count} vars: {len(variant_batch)} cons: {len(consequence_batch)} occ: {len(occurance_batch)}')
        _write_all_tables(variant_writer, conseq_writer, occurance_writer, variant_batch, consequence_batch, occurance_batch)
import pyarrow as pa
from typing import List, Dict, Any

def _process_record(case_id: int, csq_header, ped: Pedigree, record: Variant, vcf_path: str)-> tuple[list, dict, dict] | None:  
    if len(record.ALT) <= 1:
        common = process_common(record, case_id=case_id, part=0)
        picked_consequence, consequences = process_consequence(record, csq_header, common)
        occurrences = process_occurrence(record, ped, common=common)
        variant = process_variant(record, picked_consequence, common)
        return (consequences, occurrences, variant)
    else:
        logging.debug(
            f"Skipped record {record.CHROM} - {record.POS} - {record.ALT} in file {vcf_path}:"
            f" this is a multi allelic variant, mult-allelic are not supported. Please split vcf file."
        )
        return None

def _write_all_tables(
        variant_writer, conseq_writer, occurance_writer,
        variants: list[dict], consequences: list[dict], occurances: list[dict]
    ):
    _write_rows(variants, variant_schema, variant_writer)
    _write_rows(consequences, consequence_schema, conseq_writer)
    # _debug_schema_mismatch(occurances, occurance_schema)
    _write_rows(occurances, occurance_schema, occurance_writer)

def _write_rows(rows: list[dict], schema: Schema, writer):
    table = pa.Table.from_pylist(rows, schema=schema)
    writer.write_table(table)


def validate_dict_against_schema(row_dict: Dict[str, Any], schema: pa.Schema) -> List[str]:
    """
    Validate a dictionary against a schema and return error messages for mismatches.
    """
    errors = []
    
    for field in schema:
        field_name = field.name
        expected_type = field.type
        
        if field_name not in row_dict:
            if not field.nullable:
                errors.append(f"Missing required field: {field_name} (type: {expected_type})")
            continue
        
        value = row_dict[field_name]
        
        # Skip None values for nullable fields
        if value is None:
            if not field.nullable:
                errors.append(f"Null value for non-nullable field: {field_name}")
            continue
        
        # Check type compatibility
        try:
            # Try to create a scalar with the value to test type compatibility
            pa.scalar(value, type=expected_type)
        except (pa.ArrowTypeError, pa.ArrowInvalid) as e:
            errors.append(f"Field '{field_name}': Expected {expected_type}, got {type(value).__name__} with value: {value}")
    
    return errors

def _debug_schema_mismatch(rows: List[Dict[str, Any]], schema: pa.Schema):
    """
    Debug schema mismatches by checking each row and field individually.
    """
    all_errors = []
    logging.info('debug')
    for i, row in enumerate(rows):
        errors = validate_dict_against_schema(row, schema)
        if errors:
            row_errors = [f"Row {i}: {error}" for error in errors]
            all_errors.extend(row_errors)
            logging.debug(f"Errors in row {i}:")
            for error in errors:
                logging.debug(f"  {error}")
    
    return all_errors