import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
import os
import logging

def verify_parquet_file(file_path: Path, expected_schema: pa.Schema, row_count: int) -> bool:
    try:
        parquet_file = pq.ParquetFile(file_path)
        actual_schema = parquet_file.schema.to_arrow_schema()
        actual_row_count = parquet_file.metadata.num_rows
        
        if actual_schema.equals(expected_schema):
            logging.debug(f"Parquet file {file_path} exists, and schema matches")
        else:
            logging.error(f"Schema mismatch\nExpected: {expected_schema}\nActual: {actual_schema}")
            return False
        if row_count == actual_row_count:
            logging.debug(f"Row count matches ({row_count})")
            return True
        else:
            logging.debug(f"Row counts do not match. Actual: {actual_row_count} Expected: {row_count}")
            return False
            
    except Exception as e:
        logging.error(f"Error reading Parquet file: {str(e)}")
        return False