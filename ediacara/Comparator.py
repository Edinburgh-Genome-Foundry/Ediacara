# Copyright 2021 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of Ediacara.
#
# Ediacara is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Ediacara is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Ediacara. If not, see <https://www.gnu.org/licenses/>.

import itertools
import re
import statistics

import matplotlib.pyplot as plt
import numpy
import pandas
import portion as P
from weighted_levenshtein import lev

from Bio import SeqIO
import cyvcf2

import dna_features_viewer
import geneblocks


vcf_selected_info_keytable = {
    "DP": "Total read depth at the locus",
    "RO": "Reference allele observation count",
    "AO": "Alternate allele observations",
    "TYPE": "The type of allele (either snp, mnp, ins, del or complex)",
}

result_keywords = {
    "good": "PASS",
    "error": "FAIL",
    "warning": "WARNING",
    "uncertain": "LOW_COVERAGE",
}


class CustomTranslator(dna_features_viewer.BiopythonTranslator):
    """Custom translator used in the coverage plot."""

    def compute_filtered_features(self, features):
        """Display only selected features."""
        has_dnacauldron_annotation = False
        filtered_features = []
        for feature in features:
            # DNA Cauldron annotations:
            try:
                if "From " in str(feature.qualifiers.get("note", "")):
                    filtered_features += [feature]
            except:
                pass

        if len(filtered_features) == 0:  # no DNA Cauldron annotations
            # Keep up to n longest features:
            n = 10
            if len(features) > n:
                feature_lengths = []
                for feature in features:
                    feature_lengths += [len(feature.location)]
                feature_lengths.sort(reverse=True)
                max_length = feature_lengths[n]
                for feature in features:
                    if len(feature.location) > max_length:
                        filtered_features += [feature]
            else:
                filtered_features = features
        else:
            # add 4-letter overhang annotations:
            for feature in features:
                try:  # may not have a 'label'
                    if len(feature.qualifiers.get("label", "")[0]) == 4:
                        filtered_features += [feature]
                except:
                    pass

        return filtered_features


class SequencingGroup:
    """Analyse alignments of a sequencing run.


    **Parameters**

    **comparatorgroups**
    > List of `ComparatorGroup` instances.

    **name**
    > Name of the sequencing project (`str`).
    """

    def __init__(self, comparatorgroups, name="Unnamed", low_depth_cutoff=30):
        self.comparatorgroups = comparatorgroups
        self.name = name
        self.vcf_cutoff = 20  # max number of VCF entries to display in report
        self.result_keywords = (result_keywords,)
        self.low_depth_cutoff = low_depth_cutoff  # see Comparator class

    def perform_all_comparisons_in_sequencinggroup(self):
        self.number_of_reads = 0
        for comparatorgroup in self.comparatorgroups:
            comparatorgroup.perform_all_comparisons()
            self.number_of_reads += comparatorgroup.n_fastq_reads
        self.number_of_barcodes = len(self.comparatorgroups)
        self.comparisons_performed = True


class ComparatorGroup:
    """Analyse alignments to a set of references.


    **Parameters**

    **references**
    > Dictionary of name: Biopython record.

    **alignments**
    > Dictionary of pandas dataframes: {"paf": df_paf, "tsv": df_tsv}. The PAF dataframe
    was created from the output of `paftools.js sam2paf` and the TSV dataframe was
    created from the output of `samtools depth -aa`.

    **barcode**
    > The barcode number or ID (`str`).

    **assembly_paths**
    > Dictionary of sample: consensus sequence path (`dict`).

    **vcf_paths**
    > Dictionary of sample: filtered VCF filepath (`dict`).
    """

    def __init__(
        self,
        references,
        alignments,
        barcode="barcode",
        assembly_paths=None,
        vcf_paths=None,
        low_depth_cutoff=30,
    ):
        self.references = references
        self.paf = alignments["paf"]
        self.tsv = alignments["tsv"]
        self.comparators = []
        self.barcode = barcode
        self.assembly_paths = assembly_paths
        self.vcf_paths = vcf_paths
        self.low_depth_cutoff = low_depth_cutoff  # see Comparator class

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
        # ..then get only the reads aligning to reference:
        subset_paf_final = subset_paf_filtered[
            subset_paf_filtered["target_name"] == reference_name
        ]

        return subset_paf_final

    def perform_all_comparisons(self):
        """Call perform_comparison() on each Comparator.


        **Parameters**

        **assembly_paths**
        > Dictionary of construct name: consensus sequence file.
        For example: `{"Construct_1": "/path/to/assembly.fa"}`.
        """
        if self.assembly_paths is None:
            assembly_paths = {}
        else:
            assembly_paths = self.assembly_paths

        if self.vcf_paths is None:
            vcf_paths = {}
        else:
            vcf_paths = self.vcf_paths

        for comparator in self.comparators:
            try:
                assembly_path = assembly_paths[comparator.name]
            except Exception:
                assembly_path = None

            try:
                vcf_path = vcf_paths[comparator.name]
            except Exception:
                vcf_path = None

            comparator.perform_comparison(
                assembly_path=assembly_path, vcf_path=vcf_path
            )
        self.comparisons_performed = True
        self.summarise_analysis()

    def summarise_analysis(self):
        """Create variables for the PDF report summary."""
        self.number_of_constructs = len(self.comparators)
        self.result_good = 0
        self.result_warning = 0
        self.result_error = 0
        self.result_uncertain = 0

        names = []
        reference_lengths = []
        results = []
        number_of_reads_aligning = []
        median_coverages = []

        for comparator in self.comparators:
            names += [comparator.name]
            reference_lengths += [str(comparator.reference_length)]
            number_of_reads_aligning += [
                str(len(comparator.paf["query_name"].unique()))
            ]
            median_coverages += [str(int(comparator.median_yy))]

            if comparator.is_uncertain:
                self.result_uncertain += 1
                comparator.result = result_keywords["uncertain"]
            else:
                if comparator.has_errors:
                    self.result_error += 1
                    comparator.result = result_keywords["error"]
                else:
                    if comparator.has_warnings:
                        self.result_warning += 1
                        comparator.result = result_keywords["warning"]
                    else:
                        self.result_good += 1
                        comparator.result = result_keywords["good"]
            results += [comparator.result]

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
        # plt.close("all")

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

        # Mark expected lengths for interpreting the plot:
        for comparator in self.comparators:
            plt.axvline(x=len(comparator.record), color="#fd5a31")

        return fig

    @staticmethod
    def load_paf(paf_path):
        """Create a dataframe from a PAF file of alignments.


        **Parameters**

        **paf_path**
        > Path (`str`) to PAF file: `paftools.js sam2paf aln.sam > aln.paf`.
        """
        try:
            paf = pandas.read_csv(paf_path, sep="\t", header=None)
        except Exception:  # unequal number of columns, the first 12 is assumed to be OK
            # The last column is CIGAR, so we need to read lines one by one:
            rows = []
            with open(paf_path, "r") as f_input:
                paf_lines = f_input.read().splitlines()
                for paf_line in paf_lines:
                    paf_line = paf_line.split("\t")  # PAF is tab separated
                    row = paf_line[0:12] + [";".join(paf_line[12:-1])] + [paf_line[-1]]
                    rows += [row]
            paf = pandas.DataFrame(rows)
            # Set numeric columns:
            numeric_columns = [1, 2, 3, 6, 7, 8, 9, 10, 11]  # see list below
            paf[numeric_columns] = paf[numeric_columns].apply(pandas.to_numeric)

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
        self.coverage_cutoff = 0.5  # fraction of median coverage
        self.coverage_cutoff_pct = int(self.coverage_cutoff * 100)  # for PDF report
        self.reference_length = len(self.record)
        self.has_warnings = False
        self.has_errors = False
        self.geneblocks_outcome = "none"  # stores outcome, used in PDF report making
        # Set True if references are DNA Cauldron-simulated files:
        self.dnacauldron = True  # for plotting
        self.has_consensus = False  # for consensus sequence comparison
        self.is_uncertain = True  # for sequencing data quality, used in report
        self.low_depth_cutoff = 30  # for marking low-depth, value from literature

    def perform_comparison(self, assembly_path=None, vcf_path=None):
        """Plot coverage and compare reference with a consensus sequence.


        **Parameters**

        **assembly_path**
        > Optional. Path to a consensus FASTA file (`str`).
        """
        self.fig = self.plot_coverage()
        self.insert_plot = self.plot_inserts()

        if assembly_path is not None:
            self.comparison_figure = self.compare_with_assembly(
                assembly_path=assembly_path
            )
            self.has_consensus = True
        if vcf_path is not None:
            self.vcf_table = self.create_vcf_table(vcf_path=vcf_path)
        # plt.close("all")

    def create_vcf_table(self, vcf_path):
        vcf_dict = {key: [] for key in vcf_selected_info_keytable.keys()}
        vcf_dict["REF"] = []
        vcf_dict["ALT"] = []
        vcf_dict["LOC"] = []

        for variant in cyvcf2.VCF(vcf_path):
            vcf_dict["REF"] += [variant.REF]
            vcf_dict["ALT"] += [variant.ALT]
            vcf_dict["LOC"] += [variant.start]
            for key in vcf_selected_info_keytable.keys():
                vcf_dict[key] += [variant.INFO.get(key)]
        vcf_table = pandas.DataFrame(
            vcf_dict, columns=["LOC", "REF", "ALT", "TYPE", "DP", "RO", "AO"]
        )  # also set the order of the columns
        vcf_table["T"] = 1
        for index, row in vcf_table.iterrows():
            # mutations at repeats were shown to be systemic sequencing errors, and can be ignored
            # match 5 consecutive repeats:
            result = re.search(r"((\w)\2{4,})", row["REF"])
            if result is not None:
                vcf_table.loc[index, "T"] = 0
                continue

            # ignore if ALT < 50%: mutation is not real or sample is polymorphic
            if row["AO"] < row["DP"] * 0.5:
                vcf_table.loc[index, "T"] = 0

        return vcf_table

    def find_big_inserts(self, threshold=100):
        """Calculate % of reads with big inserts.


        **Parameters**

        **threshold**
        > Size of insert in bp to be considered big (`int`)."""
        read_insert_dict = {
            read_id: False for read_id in self.paf["query_name"].unique()
        }  # Initialize. Value True when read has insert.
        for index, row in self.paf.iterrows():
            read_id = row["query_name"]
            cigar = row.tail(1).item()  # CIGAR string column is the last
            for size, mutation in re.findall(r"(\d+)([A-Z]{1})", cigar):
                if int(size) >= threshold:
                    if mutation == "I":  # CIGAR code for insert
                        read_insert_dict[read_id] = True
                        break

        fraction_of_big_inserts = sum(read_insert_dict.values()) / len(
            read_insert_dict.keys()
        )
        # imprecise rounding, but ok for our purposes:
        self.pct_big_insert = round(fraction_of_big_inserts * 100, 1)  # pct multiplier

    @staticmethod
    def discretize(i):
        # see https://github.com/AlexandreDecan/portion/issues/24
        first_step = lambda s: (
            P.OPEN,
            (s.lower - 1 if s.left is P.CLOSED else s.lower),
            (s.upper + 1 if s.right is P.CLOSED else s.upper),
            P.OPEN,
        )
        second_step = lambda s: (
            P.CLOSED,
            (s.lower + 1 if s.left is P.OPEN else s.lower),
            (s.upper - 1 if s.right is P.OPEN else s.upper),
            P.CLOSED,
        )
        return i.apply(first_step).apply(second_step)

    @staticmethod
    def get_size_of_longest_interval(interval):
        if interval.empty:
            return 0

        max_length = 0
        for subinterval in interval._intervals:
            subinterval_length = subinterval.upper - subinterval.lower
            if subinterval_length > max_length:
                max_length = subinterval_length

        return max_length

    def get_unaligned_interval_sizes(self):
        grouped = self.paf.groupby("query_name")
        unaligned_interval_sizes = []
        # read_lengths = []
        for name, group in grouped:
            alignment_intervals = P.empty()
            for index, row in group.iterrows():
                interval = P.closed(row["query_start"], row["query_end"])
                alignment_intervals = alignment_intervals | interval

            read_interval = P.closed(100, row["query_length"])  # ignore first 100 bases
            discretized_alignment_intervals = self.discretize(alignment_intervals)

            unaligned_intervals = read_interval - discretized_alignment_intervals

            size = self.get_size_of_longest_interval(unaligned_intervals)
            #     read_lengths += [row["query_length"]]
            unaligned_interval_sizes += [size]

        # Get largest mode value, but handle case with no unique mode:
        occurrence = max(
            list(map(unaligned_interval_sizes.count, unaligned_interval_sizes))
        )
        self.insert_mode = max(
            list(
                set(
                    filter(
                        lambda x: unaligned_interval_sizes.count(x) == occurrence,
                        unaligned_interval_sizes,
                    )
                )
            )
        )

        self.insert_pct_above_cutoff = int(
            sum(i > 100 for i in unaligned_interval_sizes)  # 100 bp cutoff
            / len(unaligned_interval_sizes)
            * 100  # report as percent
        )
        if self.insert_pct_above_cutoff > 50:  # more than half of reads has big insert
            self.has_errors = True

        self.unaligned_interval_sizes = unaligned_interval_sizes

    def calculate_stats(self):
        """Calculate statistics for the coverage plot, used in plot_coverage()."""
        self.xx = numpy.arange(len(self.record.seq))  # for the plot x axis
        self.yy = self.tsv["depth"].to_list()  # for plotting coverage
        self.median_yy = statistics.median(self.yy)

        if self.median_yy < self.low_depth_cutoff:
            self.has_low_coverage = True
            self.has_warnings = True
        else:
            self.is_uncertain = False
            self.has_low_coverage = False

        # This section creates a list of zero coverage position to be reported
        zero_indices = [i for i, value in enumerate(self.yy) if value == 0]  # zero cov.
        G_zero = (
            list(x)
            for _, x in itertools.groupby(
                zero_indices, lambda x, c=itertools.count(): next(c) - x
            )
        )
        self.zero_coverage_positions_string = ", ".join(
            "-".join(map(str, (g[0], g[-1])[: len(g)])) for g in G_zero
        )
        if self.zero_coverage_positions_string == "":
            self.zero_coverage_positions_string = "-"  # looks better in the pdf report
        else:
            self.has_warnings = True
            self.longest_range = self.get_biggest_consecutive_range(zero_indices)
            if self.longest_range > 50:  # bp cutoff to fail sample
                self.has_errors = True

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

    def plot_coverage(self):
        """Plot the reference with the coverage and weak reads."""

        if not hasattr(self, "xx"):
            self.calculate_stats()  # also calculates yy and zz

        # Plot
        ######
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

        return fig

    def plot_inserts(self):
        """Plot the longest unaligned interval in each read (cumulative)."""

        if not hasattr(self, "unaligned_interval_sizes"):
            ## Calculate unaligned interval (insert) sizes in reads:
            self.get_unaligned_interval_sizes()

        fig, ax = plt.subplots(1, 1, figsize=(2.5, 1))
        ax.step(
            sorted(self.unaligned_interval_sizes),
            range(len(self.unaligned_interval_sizes)),
            color="red",
        )
        ax.axvline(x=100, color="grey")  # line at 100 bp for guide
        ax.set_ylabel("Reads")
        ax.set_xlabel("Interval (insert) length [bp]")

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
        """Compare the reference with a consensus sequence.

        Check length, orientation (sense or reverse complement), then use GeneBlocks to
        highlight differences between reference and the plasmid assembly.


        **Parameters**

        **assembly_path**
        > Path to a consensus FASTA file (`str`).
        """
        # Note that in this context 'assembly' means a consensus of variant calls
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
                "The consensus is %d%% (%d bp) %s than the reference."
                % (int(abs(difference_percent)), difference_bases, difference_text)
            )

            self.geneblocks_done = False
            self.is_diffblocks_reverse = None
            self.is_comparison_successful = False  # for PDF reporting
            return

        # To find out the orientation of the assembly, we compare using Levenshtein
        lev_distance = lev(str(self.assembly.seq), str(self.record.seq))
        lev_cutoff = 50  # more than this is too much to plot. Also for most plasmids,
        # this represents ~1% difference, which is too much error.
        assembly_for_diffblocks = self.assembly
        if lev_distance > lev_cutoff:
            self.perform_geneblocks = False
        else:
            self.perform_geneblocks = True

        if self.perform_geneblocks:
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

        else:  # too many differences
            self.geneblocks_done = False
            self.geneblocks_outcome = "too_many_differences"
            self.is_comparison_successful = False

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

    @staticmethod
    def get_biggest_consecutive_range(numbers):
        """Get the longest range in a list of numbers.


        **Parameters**

        **numbers**
        > A sorted list of integers (`list`).
        """
        longest = 1
        previous_number = numbers[0]
        range_length = 1
        for index, number in enumerate(numbers[1:]):
            if number == previous_number + 1:
                range_length += 1
                longest = range_length
            else:
                if range_length > longest:
                    longest = range_length
                range_length = 0
            previous_number = number

        return longest
