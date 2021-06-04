from datetime import datetime
import os

import matplotlib.pyplot as plt
import pandas

from pdf_reports import (
    dataframe_to_html,
    pug_to_html,
    write_report,
)
import pdf_reports.tools as pdf_tools

from .version import __version__

THIS_PATH = os.path.dirname(os.path.realpath(__file__))
ASSETS_PATH = os.path.join(THIS_PATH, "report_assets")
REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "simulation_report.pug")
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
    fig = comparator.plot_coverage()
    figure_data = pdf_tools.figure_data(fig, fmt="svg")
    plt.close(fig)

    comparison_figure = comparator.compare_with_assembly(assembly_path=assembly_path)
    if comparator.is_comparison_successful:
        comparison_figure_data = pdf_tools.figure_data(comparison_figure, fmt="svg")
    else:
        comparison_figure_data = None

    html = end_pug_to_html(
        REPORT_TEMPLATE,
        id=comparator.record.id,
        is_good=comparator.is_good,
        reference_length=len(comparator.record),
        figure_data=figure_data,
        low_positions=comparator.low_coverage_positions_string,
        bad_positions=comparator.high_error_positions_string,
        is_comparison_successful=comparator.is_comparison_successful,
        incorrect_length_msg=comparator.incorrect_length_msg,
        comparison_figure_data=comparison_figure_data,
        geneblocks_done=comparator.geneblocks_done,
        is_diffblocks_reverse=comparator.is_diffblocks_reverse,
    )

    # if csv_path is not None:
    #     analysis_summary.to_csv(csv_path, index=False)

    write_report(html, target, extra_stylesheets=(STYLESHEET,))
