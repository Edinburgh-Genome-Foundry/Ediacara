# Copyright 2021 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of Ediacara.
#
# Ediacara is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Ediacara is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Ediacara. If not, see <https://www.gnu.org/licenses/>.

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
