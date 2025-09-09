"""
Module for processing variant consequence annotations (VEP CSQ) from VCF records
and transforming them into a structured format compatible with Iceberg tables.

This module defines:
- A schema for consequence data.
- A helper dataclass for exon rank/total.
- Functions to parse CSQ headers, extract CSQ field values, and process consequences.

Dependencies:
- cyvcf2: for reading VCF records.
- pyiceberg: for defining the Iceberg schema.
- Common metadata and schema merging from internal modules.
"""

from cyvcf2 import Variant, VCF

from cumulus_genomic_pipeline.radiant.vcf.common import Common

CSQ_FORMAT_FIELD = "CSQ"


def get_csq_field(csq_fields, fields, field_name):
    """
    Safely retrieves a field from CSQ annotation by name.

    Args:
        csq_fields (dict): Mapping of CSQ field names to indexes.
        fields (list): Parsed CSQ annotation fields from a single transcript.
        field_name (str): Name of the CSQ field to retrieve.

    Returns:
        str or None: Value of the CSQ field or None if not found.
    """
    return fields[csq_fields[field_name]] if field_name in csq_fields else None


def process_consequence(record: Variant, csq_fields: dict[str, int], common: Common) -> tuple[dict, list[dict]]:
    """
    Processes VEP CSQ annotations from a VCF record and builds structured consequence data.

    Args:
        record (Variant): A cyvcf2 Variant object.
        csq_fields (dict[str, int]): Field name to index mapping from CSQ header.
        common (Common): Shared metadata (e.g. position, allele info).

    Returns:
        tuple:
            - dict: The primary (picked or canonical) consequence.
            - list of dict: All consequence entries for the variant.
    """
    csq = record.INFO.get(CSQ_FORMAT_FIELD, None)
    consequences = []
    pick_consequence = None
    if csq:
        csq_data = csq.split(",")
        for c in csq_data:
            fields = c.split("|")
            exon = get_csq_field(csq_fields, fields, "EXON").split("/")
            vep_impact = get_csq_field(csq_fields, fields, "IMPACT")
            hgvsg = get_csq_field(csq_fields, fields, "HGVSg")
            hgvsp = get_csq_field(csq_fields, fields, "HGVSp")
            hgvsc = get_csq_field(csq_fields, fields, "HGVSc")
            picked = get_csq_field(csq_fields, fields, "PICK") == "1"
            consequence = {
                "case_id": common.case_id,
                "locus": common.locus,
                "locus_hash": common.locus_hash,
                "chromosome": common.chromosome,
                "start": common.start,
                "end": common.end,
                "reference": common.reference,
                "alternate": common.alternate,
                "variant_class": get_csq_field(csq_fields, fields, "VARIANT_CLASS"),
                "hgvsg": hgvsg,
                "hgvsp": hgvsp,
                "hgvsc": hgvsc,
                "symbol": get_csq_field(csq_fields, fields, "SYMBOL"),
                "transcript_id": get_csq_field(csq_fields, fields, "Feature"),
                "source": get_csq_field(csq_fields, fields, "Source"),
                "biotype": get_csq_field(csq_fields, fields, "BIOTYPE"),
                "strand": get_csq_field(csq_fields, fields, "STRAND"),
                "exon": {"rank": str(exon[0]), "total": str(exon[1])} if len(exon) == 2 else None,
                "vep_impact": vep_impact,
                "consequences": get_csq_field(csq_fields, fields, "Consequence").split("&"),
                "mane_select": get_csq_field(csq_fields, fields, "ManeSelect"),
                "is_mane_select": False,
                "is_mane_plus": False,
                "is_picked": picked,
                "is_canonical": get_csq_field(csq_fields, fields, "CANONICAL") == "YES",
                "aa_change": hgvsp.split(":")[-1] if hgvsp else None,
                "dna_change": hgvsc.split(":")[-1] if hgvsp else None,
                "impact_score": IMPACT_SCORE.get(vep_impact, 0),
            }
            if picked:
                pick_consequence = consequence
            consequences.append(consequence)
    if pick_consequence is None:
        pick_consequence = next((c for c in consequences if c["is_canonical"]), None)
    return pick_consequence, consequences


def parse_csq_header(vcf: VCF):
    """
    Parses the CSQ header from a VCF and extracts field name to index mapping.

    Args:
        vcf: A cyvcf2.VCF reader object.

    Returns:
        dict[str, int]: Mapping from CSQ field names to their indexes.
    """
    info_csq = vcf.get_header_type(CSQ_FORMAT_FIELD)
    csq_meta = info_csq.get("Description", "")
    csq_meta = csq_meta.split("Format:")[-1].strip(' "')
    csq_fields = csq_meta.split("|") if csq_meta else []
    return {f: i for i, f in enumerate(csq_fields)}


IMPACT_SCORE = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
