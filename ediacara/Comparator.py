import itertools
import statistics

import matplotlib.pyplot as plt
import numpy
import pandas

from weighted_levenshtein import lev
from Bio import SeqIO

import dna_features_viewer
import geneblocks


class CustomTranslator(dna_features_viewer.BiopythonTranslator):
    """Custom translator."""

    def compute_filtered_features(self, features):
        """Display only "From " features and overhangs."""
        features = [
            feature
            for feature in features
            if (len(feature.qualifiers.get("label", "")[0]) == 4)
            or ("From " in str(feature.qualifiers.get("label", "")))
        ]
        return features


class ComparatorGroup:
    """Analyse alignments to a set of references.


    **Parameters**

    **references**
    > Dictionary of name: Biopython record.

    **alignments**
    > Dictionary of pandas dataframes: {"paf": df_paf, "tsv": df_tsv}. The PAF dataframe
    was created from the output of `paftools.js sam2paf` and the TSV dataframe was
    created from the output of `samtools depth -aa`.
    """

    def __init__(self, references, alignments):
        self.references = references
        self.paf = alignments["paf"]
        self.tsv = alignments["tsv"]
        self.comparators = []

    def create_comparator(self, reference_name):
        """Create a Comparator instance from one of the loaded references.


        **Parameters**

        **reference_name**
        > The name of the reference (`str`).
        """
        subset_paf = self.subset_paf(reference_name)
        subset_tsv = self.tsv[self.tsv["name"] == reference_name]
        alignment = {"paf": subset_paf, "tsv": subset_tsv}

        comparator = Comparator(
            {reference_name: self.references[reference_name]}, alignment
        )

        return comparator

    def add_comparator(self, reference_name):
        """Create comparator and add to the group. See create_comparator()."""
        self.comparators += [self.create_comparator(reference_name)]

    def subset_paf(self, reference_name):
        """Subset PAF dataframe for the reference.

        Select reads that align to the reference and also align *best* to the reference.


        **Parameters**

        **reference_name**
        > Name of the reference sequence to filter for (`str`).
        """
        # Work only with the values used for filtering:
        paf3 = self.paf[["query_name", "target_name", "mapping_matches"]]
        grouped = paf3.groupby("query_name")

        selected_reads = []
        for name, group in grouped:
            # not all reads align to the reference:
            if reference_name in group["target_name"].values:
                max_matches = group["mapping_matches"].max()
                # filter that the read aligns best to this reference:
                if (
                    max(
                        group[group["target_name"] == reference_name][
                            "mapping_matches"
                        ].values
                    )
                    == max_matches
                ):
                    selected_reads += [name]

        # Filter for selected reads ..
        subset_paf_filtered = self.paf[self.paf["query_name"].isin(selected_reads)]
        # .. add column for counts of multialigning reads,
        # this is for the compensation in get_weak_read_histogram_data():
        groups = (
            subset_paf_filtered.groupby("query_name")["query_name"]
            .size()
            .reset_index(name="multialign")
        )
        subset_paf_filtered_merged = subset_paf_filtered.merge(groups, how="left")
        # ..then get only the reads aligning to reference:
        subset_paf_final = subset_paf_filtered_merged[
            subset_paf_filtered_merged["target_name"] == reference_name
        ]

        return subset_paf_final

    def perform_all_comparisons(self, assembly_paths=None):
        """Call perform_comparison() on each Comparator.


        **Parameters**

        **assembly_paths**
        > Dictionary of construct name: *de novo* assembly file.
        For example: `{"Construct_1": "/path/to/assembly.fa"}`.
        """
        if assembly_paths is None:
            assembly_paths = {}
        for comparator in self.comparators:
            try:
                assembly_path = assembly_paths[comparator.name]
            except Exception:
                assembly_path = None
            comparator.perform_comparison(assembly_path=assembly_path)
        self.comparisons_performed = True
        self.summarise_analysis()

    def summarise_analysis(self):
        """Create variables for the PDF report summary."""
        self.number_of_constructs = len(self.comparators)
        self.result_good = 0
        self.result_warning = 0
        self.result_error = 0

        names = []
        reference_lengths = []
        results = []
        number_of_reads_aligning = []
        median_coverages = []

        for comparator in self.comparators:
            names += [comparator.name]
            reference_lengths += [str(comparator.reference_length)]
            number_of_reads_aligning += [str(len(comparator.paf))]
            median_coverages += [str(int(comparator.median_yy))]

            if comparator.has_errors:
                self.result_error += 1
                results += ["☒"]
            else:
                if comparator.has_warnings:
                    self.result_warning += 1
                    results += ["⚠"]
                else:
                    self.result_good += 1
                    results += ["☑"]

        self.n_fastq_reads = len(set(self.paf.query_name))
        self.fastq_plot = self.plot_fastq_histogram()

        # Summary table
        d = {
            "Name": pandas.Series(names),
            "Result": pandas.Series(results),
            "Length [bp]": pandas.Series(reference_lengths),
            "FASTQ reads": pandas.Series(number_of_reads_aligning),
            "Coverage [x]": pandas.Series(median_coverages),
        }
        self.summary_table = pandas.DataFrame(d)

    def plot_fastq_histogram(self, n_bins=50):
        """Plot a histogram of the FASTQ reads.


        **Parameters**

        **n_bins**
        > Number of bins in the histogram (`int`).
        """
        fig, ax = plt.subplots()
        plt.hist(
            self.paf["query_length"],
            bins=int(self.paf.iloc[0]["target_length"] / n_bins),
            alpha=0.8,
        )
        ax.set_ylabel("Number of reads")
        ax.set_xlabel("Read length [bp]")

        return fig

    @staticmethod
    def load_paf(paf_path):
        """Create a dataframe from a PAF file of alignments.


        **Parameters**

        **paf_path**
        > Path (`str`) to PAF file: `paftools.js sam2paf aln.sam > aln.paf`.
        """
        paf = pandas.read_csv(paf_path, sep="\t", header=None)

        # First 12 columns are defined by the format:
        columns_12 = [
            "query_name",
            "query_length",
            "query_start",
            "query_end",
            "strand",
            "target_name",
            "target_length",
            "target_start",
            "target_end",
            "mapping_matches",
            "mapping_size",
            "mapping_quality",
        ]
        # For additional columns optionally created during alignment:
        columns_13plus = [str(index) for index in range(12, len(paf.columns))]
        paf.columns = columns_12 + columns_13plus

        return paf

    @staticmethod
    def load_tsv(tsv_path):
        """Create a dataframe from a TSV file of depth counts.


        **Parameters**

        **tsv_path**
        > Path to output of `samtools depth -aa` (`str`).
        """
        tsv = pandas.read_csv(tsv_path, sep="\t", header=None)
        tsv.columns = ["name", "position", "depth"]

        return tsv


class Comparator:
    """Analyse alignments to a reference.


    **Parameters**

    **reference**
    > Dictionary of one name: one Biopython record.

    **alignment**
    > Dictionary of pandas dataframes: {"paf": df_paf, "tsv": df_tsv}. The dataframes
    contain entries only for the reference.
    """

    def __init__(self, reference, alignment):
        self.name = next(iter(reference))
        self.record = next(iter(reference.values()))
        self.paf = alignment["paf"]
        self.tsv = alignment["tsv"]
        self.coverage_cutoff = 0.15  # 15% of median coverage
        self.calculate_bad_reads = False  # for an old approach
        self.reference_length = len(self.record)
        self.has_warnings = False
        self.has_errors = False
        self.geneblocks_outcome = "none"  # stores outcome, used in PDF report making
        # Set True if references are DNA Cauldron-simulated files:
        self.dnacauldron = True  # for plotting
        self.has_de_novo = False  # for de novo assembly comparison

    def perform_comparison(self, assembly_path=None):
        """Plot coverage and compare reference with *de novo* assembly.


        **Parameters**

        **assembly_path**
        > Optional. Path to a *de novo* assembly FASTA file (`str`).
        """
        self.fig = self.plot_coverage()
        plt.close(self.fig)
        if assembly_path is not None:
            self.comparison_figure = self.compare_with_assembly(
                assembly_path=assembly_path
            )
            self.has_de_novo = True

    def calculate_stats(self):
        """Calculate statistics for the coverage plot, used in plot_coverage()."""
        self.xx = numpy.arange(len(self.record.seq))  # for the plot x axis
        self.yy = self.tsv["depth"].to_list()  # for plotting coverage
        self.median_yy = statistics.median(self.yy)

        # This section creates a list of low coverage position to be reported
        indices = [
            i
            for i, value in enumerate(self.yy)
            if value < self.median_yy * self.coverage_cutoff
        ]
        G = (
            list(x)
            for _, x in itertools.groupby(
                indices, lambda x, c=itertools.count(): next(c) - x
            )
        )
        self.low_coverage_positions_string = ", ".join(
            "-".join(map(str, (g[0], g[-1])[: len(g)])) for g in G
        )
        if self.low_coverage_positions_string == "":
            self.low_coverage_positions_string = "-"  # looks better in the pdf report
        else:
            self.has_errors = True

        # This section creates a list of bad position to be reported
        if self.calculate_bad_reads:
            if not hasattr(self, "zz"):
                self.get_weak_read_histogram_data(cutoff=0.8)
                # for plotting weak reads (less than 80% of the read maps)
            indices_errors = [
                i
                for i, value in enumerate(self.zz)
                if value > self.median_yy * self.coverage_cutoff
            ]
            H = (
                list(x)
                for _, x in itertools.groupby(
                    indices_errors, lambda x, c=itertools.count(): next(c) - x
                )
            )

            self.high_error_positions_string = ", ".join(
                "-".join(map(str, (g[0], g[-1])[: len(g)])) for g in H
            )
            if self.high_error_positions_string == "":
                self.high_error_positions_string = "-"
            else:
                self.has_errors = True

    def plot_coverage(self):
        """Plot the reference with the coverage and weak reads."""

        if not hasattr(self, "xx"):
            self.calculate_stats()  # also calculates yy and zz

        # Plot
        ######
        if self.calculate_bad_reads:
            fig, (ax1, ax2, ax3) = plt.subplots(
                3,
                1,
                figsize=(7, 4),
                sharex=True,
                gridspec_kw={"height_ratios": [4, 1, 1]},
            )
        else:
            fig, (ax1, ax2) = plt.subplots(
                2, 1, figsize=(7, 3), sharex=True, gridspec_kw={"height_ratios": [4, 1]}
            )
        if self.dnacauldron:
            graphic_record = CustomTranslator().translate_record(self.record)
        else:
            graphic_record = dna_features_viewer.BiopythonTranslator().translate_record(
                self.record
            )

        graphic_record.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=4)

        # Plot coverage
        ax2.fill_between(self.xx, self.yy, alpha=0.8)
        stdev_yy = statistics.stdev(self.yy)
        ax2.set_ylim(bottom=0)
        ax2.set_ylabel("Depth")
        ax2.set_xlabel("Base [bp]")
        ax2.axhline(
            y=self.median_yy, xmin=0, xmax=1, color="grey", linestyle="-", linewidth=0.8
        )
        ax2.axhline(
            y=self.median_yy - stdev_yy,
            xmin=0,
            xmax=1,
            color="grey",
            linestyle="--",
            linewidth=0.8,
        )
        ax2.axhline(
            y=self.median_yy + stdev_yy,
            xmin=0,
            xmax=1,
            color="grey",
            linestyle="--",
            linewidth=0.8,
        )

        # Plot low-quality reads:
        if self.calculate_bad_reads:
            ax3.fill_between(self.xx, self.zz, alpha=0.8, color="red")
            # ensure plot stays small if there are no problems:
            if max(self.zz) < self.coverage_cutoff * max(self.yy):
                ylim_top = self.coverage_cutoff * max(self.yy)
            else:
                ylim_top = max(self.zz)

            ax3.set_ylim(bottom=0, top=ylim_top)
            ax3.set_ylabel("Bad reads")
            ax3.set_xlabel("Base [bp]")
            ax3.axhline(
                y=self.median_yy * self.coverage_cutoff,
                xmin=0,
                xmax=1,
                color="grey",
                linestyle="--",
                linewidth=0.8,
            )

        return fig

    def get_weak_read_histogram_data(self, cutoff):
        """Create data for the plotting function.


        **Parameters**

        **cutoff**
        > Cutoff for proportion of bases that map (`float`). Reads below the cutoff are
        viewed as weak reads. Recommended value: 0.8 (80%).
        """
        # For each read, calculate the proportion of bases that map, then create an array
        # for each read that is below the cutoff. The sum of these arrays will give the
        # depth of aligning bad reads for each base of the reference.

        self.proportion_maps = (
            self.paf["mapping_matches"] / self.paf["query_length"].to_list()
        )
        filtered_paf = self.paf[self.proportion_maps < cutoff]
        array_length = filtered_paf.iloc[0]["target_length"]
        histogram_array = []

        for index, row in filtered_paf.iterrows():
            read_array = [0] * array_length
            start = row["target_start"]
            end = row["target_end"]
            value = 1 / row["multialign"]
            values = [value] * (end - start)
            read_array[start:end] = values
            histogram_array.append(read_array)

        histogram_numpy_array = numpy.array(histogram_array)

        self.zz = numpy.sum(histogram_numpy_array, axis=0)

    def compare_with_assembly(self, assembly_path):
        """Compare the reference with a Canu assembly file.

        Check length, orientation (sense or reverse complement), then use GeneBlocks to
        highlight differences between reference and the plasmid assembly.


        **Parameters**

        **assembly_path**
        > Path to Canu assembly FASTA file (`str`).
        """
        # Note that in this context 'assembly' means a de novo genome assembly created
        # from sequencing reads.
        self.assembly = SeqIO.read(handle=assembly_path, format="fasta")
        length_ratio = len(self.assembly) / len(self.record)
        if 0.95 < length_ratio < 1.05:  # length within reason
            self.is_length_ok = True  # can compare with expected sequence
            self.incorrect_length_msg = None
        else:
            self.geneblocks_outcome = "incorrect_length"
            self.is_length_ok = False  # do not compare sequences
            self.has_errors = True

            difference_bases = len(self.assembly) - len(self.record)
            difference_percent = (length_ratio - 1) * 100  # get as percent
            if difference_percent > 0:
                difference_text = "longer"
            else:
                difference_text = "shorter"
            self.incorrect_length_msg = (
                "The <i>de novo</i> assembly is %d%% (%d bp) %s than the reference."
                % (int(abs(difference_percent)), difference_bases, difference_text)
            )

            self.geneblocks_done = False
            self.is_diffblocks_reverse = None
            self.is_comparison_successful = False  # for PDF reporting
            return

        # To find out the orientation of the assembly, we compare using Levenshtein
        lev_distance = lev(str(self.assembly.seq), str(self.record.seq))
        lev_rc_distance = lev(
            str(self.assembly.seq.reverse_complement()), str(self.record.seq)
        )
        if lev_distance < lev_rc_distance:
            self.is_assembly_reverse_complement = False
            assembly_for_diffblocks = self.assembly
        else:
            self.is_assembly_reverse_complement = True  # for geneblocks
            assembly_for_diffblocks = self.assembly.reverse_complement()

        try:
            diff_blocks = geneblocks.DiffBlocks.from_sequences(
                assembly_for_diffblocks, self.record
            )
        except KeyError:
            self.geneblocks_done = False
            self.geneblocks_outcome = "geneblocks_error"
        else:  # diffblocks had no error
            self.geneblocks_done = True
            ax1, ax2 = diff_blocks.plot(figure_width=5)
            self.is_diffblocks_reverse = False
            self.is_comparison_successful = True
            self.geneblocks_outcome = "all_good"

            return ax1

        # We try again as due to a bug in geneblocks, the reverse order may work:
        if self.geneblocks_outcome == "geneblocks_error":
            try:
                diff_blocks = geneblocks.DiffBlocks.from_sequences(
                    self.record, assembly_for_diffblocks
                )
            except KeyError:
                return

            else:
                self.geneblocks_done = True  # overwrite
                ax1, ax2 = diff_blocks.plot(figure_width=7)
                self.is_diffblocks_reverse = True
                self.is_comparison_successful = True
                self.geneblocks_outcome = "swapped_diffblocks"

                return ax1

    def filter_fastq(self, fastq_path, target=None):
        """Filter original FASTQ file for reads that best map to the reference.

        This function is used to obtain a FASTQ file that can be used for Canu assembly.


        **Parameters**

        **fastq_path**
        > Path to FASTQ file (`str`).

        **target**
        > Path to save the filtered FASTQ file (`str`).
        """
        if target is None:
            target = fastq_path + "_filtered.fastq"
        read_names = self.paf["query_name"].to_list()
        input_seq_iterator = SeqIO.parse(fastq_path, "fastq")
        seq_iterator = (
            record for record in input_seq_iterator if record.name in read_names
        )
        SeqIO.write(seq_iterator, target, "fastq")
