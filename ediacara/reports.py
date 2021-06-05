from datetime import datetime
import os

import matplotlib.pyplot as plt
import pandas

from pdf_reports import (
    add_css_class,
    dataframe_to_html,
    pug_to_html,
    style_table_rows,
    write_report,
)
import pdf_reports.tools as pdf_tools

from .version import __version__

THIS_PATH = os.path.dirname(os.path.realpath(__file__))
ASSETS_PATH = os.path.join(THIS_PATH, "report_assets")
REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "simulation_report.pug")
GROUP_REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "group_simulation_report.pug")
STYLESHEET = os.path.join(ASSETS_PATH, "report_style.css")


def end_pug_to_html(template, **context):
    now = datetime.now().strftime("%Y-%m-%d")
    defaults = {
        "sidebar_text": "Generated on %s by EGF's alignment comparator (version %s)"
        % (now, __version__),
        "end_logo_url": os.path.join(ASSETS_PATH, "imgs", "logo.png"),
    }
    for k in defaults:
        if k not in context:
            context[k] = defaults[k]
    return pug_to_html(template, **context)


def write_comparatorgroup_report(target, comparatorgroup):
    """Write an alignment report with a PDF summary.


    **Parameters**

    **target**
    > Path for PDF file.

    **comparatorgroup**
    > ComparatorGroup instance.
    """
    if not comparatorgroup.comparisons_performed:
        return "Run perform_all_comparisons()!"

    comparatorgroup.report_table = dataframe_to_html(
        comparatorgroup.summary_table, extra_classes=("definition",)
    )

    def tr_modifier(tr):
        tds = list(tr.find_all("td"))
        if len(tds) == 0:
            return
        result = tds[1]  # second element of list is the result symbol
        if result.text == "☑":
            add_css_class(tr, "positive")
        elif result.text == "☒":
            add_css_class(tr, "negative")

    # This colours the summary table:
    comparatorgroup.report_table = style_table_rows(
        comparatorgroup.report_table, tr_modifier
    )

    # Histogram of reads in the report summary
    comparatorgroup.fastq_figure_data = pdf_tools.figure_data(
        comparatorgroup.fastq_plot, fmt="svg"
    )

    for comparator in comparatorgroup.comparators:
        comparator.figure_data = pdf_tools.figure_data(comparator.fig, fmt="svg")

        if hasattr(comparator, "is_comparison_successful"):
            if comparator.is_comparison_successful:
                comparator.has_comparison = True
                height = comparator.comparison_figure.figure.get_size_inches()[1]
                if height > 10:
                    height = 10  # to fit on one page
                comparator.comparison_figure_data = pdf_tools.figure_data(
                    comparator.comparison_figure, fmt="svg", size=[7, height]
                )
            else:
                comparator.comparison_figure_data = None

    html = end_pug_to_html(GROUP_REPORT_TEMPLATE, comparatorgroup=comparatorgroup)
    write_report(html, target, extra_stylesheets=(STYLESHEET,))


def write_pdf_report(
    target, comparator, assembly_path, csv_path=None,
):
    """Write an alignment report with a PDF summary.


    **Parameters**

    **target**
    > Path for PDF file.

    **comparator**
    > Comparator instance.

    **assembly_path**
    > Assembly to compare with. See `Comparator.compare_with_assembly()`.

    **csv_path**
    > Path to CSV output of results (`str`). Default no CSV.
    """
    comparator.perform_comparison(assembly_path)
    comparator.figure_data = pdf_tools.figure_data(comparator.fig, fmt="svg")

    if hasattr(comparator, "is_comparison_successful"):
        if comparator.is_comparison_successful:
            comparator.has_comparison = True
            comparator.comparison_figure_data = pdf_tools.figure_data(
                comparator.comparison_figure, fmt="svg"
            )
        else:
            comparator.comparison_figure_data = None

    html = end_pug_to_html(REPORT_TEMPLATE, comparator=comparator,)

    # if csv_path is not None:
    #     analysis_summary.to_csv(csv_path, index=False)

    write_report(html, target, extra_stylesheets=(STYLESHEET,))
