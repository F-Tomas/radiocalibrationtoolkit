#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import skrf as rf
from pathlib import Path
from radiocalibrationtoolkit import *

MHz = 1e6


def nparray2string(arr):
    string = ""
    for i in arr:
        # string+='{:.2f} '.format(i)
        string += "{} ".format(round(i, 2))
    return string


def get_hw_profile_xml_string(
    responseId,
    freq_string,
    mag_string,
    phase_string,
    mag_unit="dB",
    phase_unit="degree",
    type_="Response",
):
    freq_string = nparray2string(freq_string)
    mag_string = nparray2string(mag_string)

    header = "<!-- {} -->".format(ifilename)

    freqUnit = "MHz"

    type__ = type_[0] + type_

    if phase_string is not None:
        phase_string = nparray2string(phase_string)
        phase_XML_string = (
            '      <Phase unit ="' + phase_unit + '">\n'
            "        " + phase_string + "      </Phase>\n"
        )
    else:
        phase_XML_string = ""

    xml_string = (
        header
        + "\n"
        + "  <"
        + type__
        + ' id="'
        + responseId
        + '"> \n'
        + "    <"
        + type_
        + "> \n "
        + '      <x unit="'
        + freqUnit
        + '"> \n'
        + "        "
        + freq_string
        + "      </x>\n"
        + "      <"
        + mag_unit
        + ">\n"
        + "        "
        + mag_string
        + "      </"
        + mag_unit
        + ">\n"
        + phase_XML_string
        + "    </"
        + type_
        + ">\n"
        + "  </"
        + type__
        + ">"
        + "\n\n\n"
    )
    return xml_string


def calculate_impedance(c):
    return args["reference_impedance"] * (1 + c) / (1 - c)


def s2impedance(sparameter_file_path, frequency_MHz_new=None):
    ntwk = rf.Network(sparameter_file_path)
    if frequency_MHz_new != None:
        ntwk = ntwk.interpolate(frequency_MHz_new, fill_value="extrapolate")
    return ntwk.f / MHz, np.real(
        np.asarray(list(map(calculate_impedance, ntwk.s11.s.flatten())))
    )


def s2dB(sparameter_file_path, unwrapped_phase=True, frequency_MHz_new=None):
    ntwk = rf.Network(sparameter_file_path)
    if frequency_MHz_new != None:
        ntwk = ntwk.interpolate(frequency_MHz_new, fill_value="extrapolate")

    frequency_MHz = ntwk.f / MHz
    magnitude_dB = ntwk.s21.s_db.flatten()
    phase_degs_wrapped = ntwk.s21.s_deg.flatten()

    # unwrape phase if wraped
    phase_degs_unwrapped = np.unwrap(phase_degs_wrapped, 180)

    if unwrapped_phase:
        phase_degs = phase_degs_unwrapped
    else:
        phase_degs = phase_degs_wrapped

    magnitude_dB = np.around(magnitude_dB, args["round"])
    phase_degs = np.around(phase_degs, args["round"])

    return frequency_MHz, magnitude_dB, phase_degs


def main():
    if args["frange"] != None:
        frequency_MHz_new = rf.Frequency(
            args["frange"][0],
            args["frange"][1],
            args["frange"][1] - args["frange"][0] + 1,
            "MHz",
        )
    else:
        frequency_MHz_new = None

    frequency_MHz, values_dB, values_phase_degs = s2dB(
        args["inputfile"], frequency_MHz_new=frequency_MHz_new
    )
    amp_xml_string = get_hw_profile_xml_string(
        label,
        frequency_MHz,
        values_dB,
        values_phase_degs,
        mag_unit="dB",
        phase_unit="degree",
    )

    frequency_MHz, impedance = s2impedance(
        args["inputfile"], frequency_MHz_new=frequency_MHz_new
    )
    imp_string = get_hw_profile_xml_string(
        label,
        frequency_MHz,
        impedance,
        None,
        mag_unit="y",
        phase_unit="degree",
        type_="Impedance",
    )

    abs_output_path = os.path.join(args["outputpath"], args["savefilename"])
    closing_tags = "</HardwareProfileList>"
    # if appending, delete closing tags and do not use the opening tags
    if args["mode"] == "a":
        opening_tags = ""
        with open(abs_output_path, "r+") as f:
            file_contents = f.read()
            file_contents = file_contents.replace(closing_tags, "")
            f.seek(0)
            f.write(file_contents)
            f.truncate()
    else:
        opening_tags = "<HardwareProfileList>\n\n"

    with open(abs_output_path, args["mode"]) as f:
        f.write(opening_tags + amp_xml_string + imp_string + closing_tags)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap._action_groups.pop()
    required = ap.add_argument_group("required arguments")
    optional = ap.add_argument_group("optional arguments")

    required.add_argument(
        "-i",
        "--inputfile",
        nargs="?",
        required=True,
        help="Path to the S-parameter file.",
    )
    required.add_argument(
        "-m",
        "--mode",
        nargs="?",
        required=True,
        choices=["w", "a"],
        help="Append or Write to an existing file.",
    )
    optional.add_argument(
        "-l",
        "--label",
        nargs="?",
        default="",
        help="Label for the hardware. If not provided, stem from the file name will be used.",
    ),
    optional.add_argument(
        "-f",
        "--frange",
        nargs="+",
        default=None,
        type=int,
        help="Use skrf to interpolate the s-parameter file values to a frequency range in MHz, 'from' 'to'.\
              Must be an integer!. Spacing is '1'",
    ),
    optional.add_argument(
        "-o",
        "--outputpath",
        nargs="?",
        default="./antenna_setup_files/",
        help="Output folder. If it does not exists, it will be created.",
    )
    optional.add_argument(
        "-s",
        "--savefilename",
        nargs="?",
        default="HardwareProfileList.xml",
        help="Saved power dataset filename.",
    )
    optional.add_argument(
        "-r",
        "--round",
        nargs=1,
        default=3,
        type=int,
        help="To how many places the hw profiles will be rounded..",
    )
    optional.add_argument(
        "-z",
        "--reference_impedance",
        nargs=1,
        default=50,
        type=float,
        help="Reference impedance, default is 50 Ohms.",
    )

    args = vars(ap.parse_args())

    ifilename = Path(args["inputfile"]).stem
    if args["label"] == "":
        label = ifilename
    else:
        label = args["label"]

    mkdir(args["outputpath"])

    sys.exit(main())
