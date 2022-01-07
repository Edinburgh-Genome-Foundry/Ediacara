from Bio import SeqIO

from .Comparator import ComparatorGroup


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
