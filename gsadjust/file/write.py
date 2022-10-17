"""
file/write.py
===============

Data export code for GSadjust.
--------------------------------------------------------------------------------

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""
import csv
import logging
import os
import time

from ..data.analysis import compute_gravity_change


SEPARATOR = " | "


def line_generator(lines):
    for line in lines:
        yield line + "\n"


def export_metadata(obsTreeModel, data_path):
    """
    Write metadata text to file. Useful for USGS data releases.

    The filename is generated automatically with a timestamp.

    Parameters
    ----------
    obsTreeModel : MainProg data object
    data_path : str
        Path to write output

    Returns
    -------
    filename or False
    """
    filename = os.path.join(
        data_path, "GSadjust_MetadataText_" + time.strftime("%Y%m%d-%H%M") + ".txt"
    )

    output_format = "table"

    export_fn = {
        "table": _export_metadata_table,
        "text": _export_metadata_text,
    }

    fn = export_fn.get(output_format)
    lines = fn(obsTreeModel)

    if not lines:
        return False

    with open(filename, "w") as fid:
        fid.write(
            "Attribute accuracy is evaluated from the least-squares network adjustment"
            " results.\n"
        )
        fid.writelines(lines)

    return filename


def _export_metadata_table(obsTreeModel):
    table = [
        SEPARATOR.join(
            [
                "Survey",
                "Max. delta residual",
                "Max. datum residual",
                "Mean SD",
                "Deltas",
                "Deltas not used",
                "Datums",
                "Datums not used",
            ]
        )
    ]
    sf_header_written = False
    for survey in obsTreeModel.checked_surveys():
        if survey.adjustment.adjustmentresults.n_datums > 0:
            table.append(
                f"{survey.name} "
                f"{survey.adjustment.adjustmentresults.max_dg_residual:>5.1f} "
                f"{survey.adjustment.adjustmentresults.max_datum_residual:>5.1f} "
                f"{survey.adjustment.adjustmentresults.avg_stdev:>5.1f} "
                f"{survey.adjustment.adjustmentresults.n_deltas:>4} "
                f"{survey.adjustment.adjustmentresults.n_deltas_notused:>3} "
                f"{survey.adjustment.adjustmentresults.n_datums:>3} "
                f"{survey.adjustment.adjustmentresults.n_datums_notused:>3}"
            )
    for survey in obsTreeModel.checked_surveys():
        if survey.adjustment.adjustmentresults.n_datums > 0:
            if len(survey.adjustment.adjustmentresults.cal_dic) > 0:
                if not sf_header_written:
                    table.append("Relative gravimeter scale factor(s)")
                    table.append(
                        "Survey | Meter | Scale factor | Scale factor S.D. (0 ="
                        " specified S.F.)"
                    )
                    sf_header_written = True
                for k, v in survey.adjustment.adjustmentresults.cal_dic.items():
                    table.append(
                        "{} {:>6} {:>10.6f} {:>10.6f}".format(
                            survey.name, k, v[0], v[1]
                        )
                    )
            elif survey.adjustment.adjustmentoptions.specify_cal_coeff:
                if not sf_header_written:
                    table.append("Relative gravimeter scale factor(s)")
                    table.append("Survey | Meter | Scale factor | Scale factor S.D.")
                    sf_header_written = True
                for (
                    k,
                    v,
                ) in survey.adjustment.adjustmentoptions.meter_cal_dict.items():
                    table.append("{} {:>6} {:>10.6f} 0".format(survey.name, k, v))
    return line_generator(table) if table else False


def _export_metadata_text(obsTreeModel):
    lines = []
    for survey in obsTreeModel.checked_surveys():
        if (
            survey.adjustment.adjustmentresults.n_datums > 0
        ):  # check that there are results
            lines.append(
                f"For the {survey.name} survey, the minimum and maximum"
                " gravity-difference residuals were"
                f" {survey.adjustment.adjustmentresults.min_dg_residual:0.1f} and"
                f" {survey.adjustment.adjustmentresults.max_dg_residual:0.1f} microGal,"
                " respectively. The minimum and maximum datum (absolute-gravity"
                " station) residuals were"
                f" {survey.adjustment.adjustmentresults.min_datum_residual:0.1f} and"
                f" {survey.adjustment.adjustmentresults.max_datum_residual:0.1f} microGal,"
                " respectively. The average standard deviation of the adjusted gravity"
                " values at each station (derived from the network adjustment) was"
                f" {survey.adjustment.adjustmentresults.avg_stdev:0.1f} microGal. "
            )
            # TODO: account for instance of 1 outlier ('1 was removed')
            datum_was_or_were, delta_was_or_were = "were", "were"
            outlier_or_outliers = "outliers"
            if survey.adjustment.adjustmentresults.n_datums == 1:
                datum_was_or_were = "was"
            if survey.adjustment.adjustmentresults.n_deltas_notused == 1:
                delta_was_or_were = "was"
                outlier_or_outliers = "an outlier"
            lines.append(
                "{} out of {} possible gravity differences were used in the adjustment"
                " ({} {} removed ".format(
                    survey.adjustment.adjustmentresults.n_deltas,
                    survey.adjustment.adjustmentresults.n_deltas_notused
                    + survey.adjustment.adjustmentresults.n_deltas,
                    survey.adjustment.adjustmentresults.n_deltas_notused,
                    delta_was_or_were,
                )
            )
            lines.append(
                "as {}). {} out of {} possible datum observations {} used. ".format(
                    outlier_or_outliers,
                    survey.adjustment.adjustmentresults.n_datums,
                    survey.adjustment.adjustmentresults.n_datums_notused
                    + survey.adjustment.adjustmentresults.n_datums,
                    datum_was_or_were,
                )
            )
            logging.info("Metadata text written to file")

    return line_generator(lines) if lines else False


def export_summary(obsTreeModel, data_path):
    """
    Write summary of procesing to text file. Can be used to reproduce results.

    The filename is generated automatically with a timestamp.

    Parameters
    ----------
    obsTreeModel : MainProg data object
    data_path : str
        Path to write output

    Returns
    -------
    filename or False
    """
    fn = os.path.join(
        data_path, "GSadjust_Summary_" + time.strftime("%Y%m%d-%H%M") + ".txt"
    )
    # Write header info
    with open(fn, "w") as fid:
        fid.write(
            "# GSadjust processing summary, {}\n".format(
                time.strftime("%Y-%m-%d %H:%M")
            )
        )
        fid.write("# Station data\n")
        for survey in obsTreeModel.surveys():
            for loop in survey.loops():
                fid.write(
                    "# Survey {}, Loop {}, Source file: {}\n".format(
                        survey.name, loop.name, loop.source
                    )
                )
                for iii in range(loop.rowCount()):
                    station = loop.child(iii)
                    fid.write("# " + station.summary_str)
                    fid.write(
                        "# Checked | Station | Raw gravity | ET correction | Corr."
                        " gravity | Std. dev.\n"
                    )
                    for sample_str in station.iter_samples():
                        fid.write(sample_str)
        fid.write("# Loop data\n")
        for survey in obsTreeModel.surveys():
            for ii in range(survey.rowCount()):
                loop = survey.child(ii)
                fid.write(
                    "# Survey {}, Loop {}, Source file: {}\n# ".format(
                        survey.name, loop.name, loop.source
                    )
                )
                fid.write(str(loop))
                fid.write(
                    "# Checked | Station1 | Station2 | Date | Time (UTC) | delta-g |"
                    " Std. dev. | Drift correction\n"
                )
                for delta in loop.deltas:
                    fid.write(
                        "{} {} {} {} {:.2f} {:.2f} {}\n".format(
                            int(delta.checked / 2),
                            delta.sta1,
                            delta.sta2,
                            delta.time_string(),
                            delta.dg,
                            delta.sd,
                            delta.driftcorr,
                        )
                    )
        fid.write("# Adjustment data\n")
        for survey in obsTreeModel.surveys():
            if survey.checkState() == 0:
                fid.write(
                    "Survey {} not included in results (survey was unchecked)\n".format(
                        survey.name
                    )
                )
            else:
                fid.write("# Adjustment options, survey: {}\n".format(survey.name))
                fid.write(str(survey.adjustment.adjustmentoptions))
                fid.write(
                    "# Deltas in the adjustment, survey: {}\n".format(survey.name)
                )
                fid.write(
                    "# Checked | Station1 | Station2 | Date | Time (UTC) | delta-g |"
                    " Std. dev. | Std. dev. for adj. | Residual\n"
                )
                for delta in survey.adjustment.deltas:
                    fid.write("{}\n".format(delta))
                fid.write(
                    "# Datums in the adjustment, survey: {}\n".format(survey.name)
                )
                fid.write(
                    "# Checked | Station | Date | Gravity | Std. dev. | Residual\n"
                )
                for datum in survey.adjustment.datums:
                    fid.write("{}\n".format(datum))
                fid.write("# Adjustment results, survey: {}\n".format(survey.name))

                lines = line_generator(survey.adjustment.results_string())
                fid.writelines(lines)

                fid.write("# Adjusted station values, survey: {}\n".format(survey.name))
                fid.write("# Station | Gravity | Std. dev.\n")
                for adj_sta in survey.results:
                    fid.write("{}\n".format(str(adj_sta)))
    return fn


def export_data(obstreemodel, data_path):
    """
    Export gravity change table to csv file

    The filename is generated automatically with a timestamp.

    Parameters
    ----------
    obstreemodel : MainProg data object
    data_path : str
        Path to write output

    Returns
    -------
    filename or False

    """
    fn = os.path.join(
        data_path, "GSadjust_TabularData_" + time.strftime("%Y%m%d-%H%M") + ".csv"
    )
    table = compute_gravity_change(obstreemodel, table_type="full")

    with open(fn, "w", newline="\n") as fid:
        wr = csv.writer(fid)
        wr.writerow(table[0])
        for row in table[1]:
            wr.writerow(row)

    return fn
