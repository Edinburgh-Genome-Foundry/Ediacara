import os

from Bio import SeqIO

import ediacara as edi

data_dir = os.path.join("tests", "data")


# def test_comparator():
#     edi.Comparator


def test_comparatorgroup():
    # Also tests Comparator
    sample = "EGF2_2"
    vcf = os.path.join(data_dir, "barcode01_EGF2_2_filtered.vcf")

    reference_gb = os.path.join(data_dir, "EGF2_2.gb")
    record = SeqIO.read(reference_gb, "genbank")
    record.id = sample
    references = {record.id: record}

    tsv_file = os.path.join(data_dir, "barcode01_EGF2_2_sorted.sam_counts.tsv")
    paf_path = os.path.join(data_dir, "barcode01_EGF2_2.paf")

    tsv = edi.ComparatorGroup.load_tsv(tsv_file)
    paf = edi.ComparatorGroup.load_paf(paf_path)

    assembly_paths = {sample: os.path.join(data_dir, "EGF2_2.fa")}
    vcf_paths = {sample: vcf}

    comparator_group = edi.ComparatorGroup(
        references=references,
        alignments={"paf": paf, "tsv": tsv},
        barcode="barcode01",
        assembly_paths=assembly_paths,
        vcf_paths=vcf_paths,
    )

    list_of_constructs = [sample]
    for element in list_of_constructs:
        comparator_group.add_comparator(element)
