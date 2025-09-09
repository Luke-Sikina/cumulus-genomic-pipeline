# Cumulus Genomic Pipeline

This is a primitive pipeline that processes VCFs into parquet files. It uses the RADIANT portal parquet schema and ETL logic.

## Requirements

- python3
- poetry
- DuckDB (for viewing outputs)

## Usage

Build and run:
```shell
poetry install
poetry run python src/cumulus_genomic_pipeline/main.py -v -i tests/data/4klines.variants.CEPH-1463.snv.vep.vcf.gz -o out/
```

View results:
```shell
duckdb
#D select * from 'out/variants.parquet';
```

## Development

Build: `poetry install`

Run tests: `poetry run pytest`