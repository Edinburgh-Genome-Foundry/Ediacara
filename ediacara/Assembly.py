from Bio import SeqFeature
from Bio import SeqIO

import dna_features_viewer

from .Comparator import ComparatorGroup


class AssemblyTranslator(dna_features_viewer.BiopythonTranslator):
    """Custom translator to display only the aligned parts."""

    def compute_feature_color(self, feature):
        ediacara_qualifier = "ediacara"
        if feature.qualifiers[ediacara_qualifier] == "wrong_part":
            return "#f5685e"
        elif feature.qualifiers[ediacara_qualifier] == "correct_part":
            return "#79d300"
        else:
            return "#7245dc"  # default dna_features_viewer colour


class Assembly:
    """Compare an assembly sequence to parts and simulated reference.


    **Parameters**

    **assembly_path**
    > Path to the *de novo* assembled FASTA sequence file (`str`).

    **reference**
    > Path to reference Genbank file (`str`).
    
    **alignment**
    > Path (`str`) to minimap2 alignment PAF file, created with the `-cx asm5` options.
    In this, the part and reference sequences are aligned against the de novo sequence.
    """

    def __init__(self, assembly_path, reference_path, alignment_path):
        self.assembly = SeqIO.read(handle=assembly_path, format="fasta")
        self.reference = SeqIO.read(handle=reference_path, format="genbank")
        self.paf = ComparatorGroup.load_paf(alignment_path)
        self.paf.columns = self.paf.columns[:-1].to_list() + ["CIGAR"]  # -1 is last

    def interpret_alignment(self):
        for index, row in self.paf.iterrows():
            # May be useful to exclude reference:
            # if row["query_name"] == self.reference.id:
            #     continue
            location = SeqFeature.FeatureLocation(
                row["target_start"], row["target_end"]
            )
            feature = SeqFeature.SeqFeature(
                location=location,
                type="misc_feature",
                id=row["query_name"],
                qualifiers={"label": row["query_name"], "ediacara": "correct_part"},
            )

            self.assembly.features.append(feature)
