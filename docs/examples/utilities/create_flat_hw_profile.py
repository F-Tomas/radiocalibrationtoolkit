#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import sys

def nparray2string(arr):
    string = ""
    for i in arr:
        # string+='{:.2f} '.format(i)
        string += "{} ".format(round(i, 2))
    return string


def get_hw_profile_xml_string(
    responseId,
    freqString,
    magString,
    phaseString,
    magUnit="dB",
    phaseUnit="degree",
    type_="Response",
):
    freqString = nparray2string(freqString)
    magString = nparray2string(magString)

    header = "<!-- {} -->".format(responseId)

    freqUnit = "MHz"
    magUnit = "dB"
    phaseUnit = "degree"

    type__ = type_[0] + type_

    if phaseString is not None:
        phase_XML_string = (
            '      <Phase unit ="' + phaseUnit + '">\n'
            "        " + phaseString + "      </Phase>\n"
        )
    else:
        phase_XML_string = ""

    xml_string = (
        header
        + "\n\n\n"
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
        + freqString
        + "      </x>\n"
        + "      <"
        + magUnit
        + ">\n"
        + "        "
        + magString
        + "      </"
        + magUnit
        + ">\n"
        + phase_XML_string
        + "    </"
        + type_
        + ">\n"
        + "  </"
        + type__
        + ">"
    )
    return xml_string


def main():
    flat_hw_0 = -6.35
    flat_hw_1 = 13.12
    flat_hw_2 = -1.05
    flat_hw_3 = 17.64
    flat_impedance = 400

    frequency = np.asarray(list(range(1, 126)))
    custom_hw_0 = np.ones(len(frequency)) * flat_hw_0
    custom_hw_1 = np.ones(len(frequency)) * flat_hw_1
    custom_hw_2 = np.ones(len(frequency)) * flat_hw_2
    custom_hw_3 = np.ones(len(frequency)) * flat_hw_3
    # Supression until 30 MHz and from 80 MHz
    custom_hw_3[(frequency < 30) | (frequency > 80)] = -20

    hw_xml_string_0 = get_hw_profile_xml_string(
        "impedance_matching_EW",
        frequency,
        custom_hw_0,
        None,
        magUnit="dB",
        phaseUnit="degree",
    )
    hw_xml_string_1 = get_hw_profile_xml_string(
        "LNA",
        frequency,
        custom_hw_1,
        None,
        magUnit="dB",
        phaseUnit="degree",
    )
    hw_xml_string_2 = get_hw_profile_xml_string(
        "cable_fromLNA2digitizer",
        frequency,
        custom_hw_2,
        None,
        magUnit="dB",
        phaseUnit="degree",
    )
    hw_xml_string_3 = get_hw_profile_xml_string(
        "digitizer",
        frequency,
        custom_hw_3,
        None,
        magUnit="dB",
        phaseUnit="degree",
    )

    imp_string = get_hw_profile_xml_string(
        "antenna_EW",
        frequency,
        np.ones(frequency.size) * flat_impedance,
        None,
        magUnit="y",
        type_="Impedance",
    )

    hw_xml_string = (
        hw_xml_string_0
        + "\n\n\n"
        + hw_xml_string_1
        + "\n\n\n"
        + hw_xml_string_2
        + "\n\n\n"
        + hw_xml_string_3
        + "\n\n\n"
        + imp_string
        + "\n\n\n"
    )

    closing_tags = "</HardwareProfileList>"
    opening_tags = "<HardwareProfileList>\n\n"

    with open("./antenna_setup_files/HardwareProfileList_flat.xml", "w") as file:
        file.write(opening_tags+hw_xml_string+closing_tags)


if __name__ == "__main__":
    sys.exit(main())
