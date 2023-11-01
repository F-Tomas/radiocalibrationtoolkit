#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import argparse


def nparray2string(arr):
    string = ""
    for i in arr:
        # string+='{:.2f} '.format(i)
        string += "{} ".format(round(i, 2))
    return string


def main():
    thetaList, phiList = np.meshgrid(
        np.asarray(list(range(0, 91, 5))), np.asarray(list(range(0, 361, 5)))
    )
    thetaList = thetaList.flatten()
    phiList = phiList.flatten()

    freqList = np.asarray(list(range(0, 126, 1)))
    vel_shifted_mean = np.ones(freqList.size)

    vel_theta = {}
    vel_phi = {}
    vel_ephi_phase = {}
    vel_etheta_phase = {}
    for f in freqList:
        vel_theta[f] = np.ones(thetaList.size)
        vel_phi[f] = np.ones(phiList.size)
        vel_ephi_phase[f] = np.zeros(phiList.size)
        vel_etheta_phase[f] = np.zeros(thetaList.size)

    # antennaPatternXMLfile = "../antenna_setup_files/Isotropic_antenna_pattern.xml"

    with open(antennaPatternXMLfile, "w") as xml_file:
        xml_file.write('<frequency unit="MHz">  <!-- Simulated Frequencies -->\n')
        xml_file.write(nparray2string(freqList) + "\n")
        xml_file.write(
            '</frequency>\n<theta unit="deg">  <!-- zenith angle theta -->\n'
        )
        xml_file.write(nparray2string(thetaList) + "\n")
        xml_file.write('</theta>\n<phi unit="deg">  <!-- azimuth angle phi -->\n')
        xml_file.write(nparray2string(phiList) + "\n")
        xml_file.write(
            '</phi>\n<MeanTransfer unit="m">  <!-- Mean Transfer of the antenna to electric field - not directional -->\n'
        )
        xml_file.write(nparray2string(vel_shifted_mean) + "\n")
        xml_file.write("</MeanTransfer>\n")
        #
        xml_file.write(
            "<!-- EAHTheta_amp: Projection of the effective antenna height in e_theta direction -->\n"
        )
        for f in freqList:
            xml_file.write('<EAHTheta_amp unit="m" idfreq="%.1f">' % f)
            xml_file.write(nparray2string(vel_theta[f]) + "</EAHTheta_amp>\n")
        xml_file.write("<!-- end of EAHTheta_amp -->\n")
        #
        xml_file.write(
            "<!-- EAHTheta_phase: Phase of the effective antenna height in e_theta direction -->\n"
        )
        for f in freqList:
            xml_file.write('<EAHTheta_phase unit="deg" idfreq="%.1f">' % f)
            xml_file.write(nparray2string(vel_etheta_phase[f]) + "</EAHTheta_phase>\n")
        xml_file.write("<!-- end of EAHTheta_phase -->\n")
        #
        xml_file.write(
            "<!-- EAHPhi_amp: Projection of the effective antenna height in e_phi direction -->\n"
        )
        for f in freqList:
            xml_file.write('<EAHPhi_amp unit="m" idfreq="%.1f">' % f)
            xml_file.write(nparray2string(vel_phi[f]) + "</EAHPhi_amp>\n")
        xml_file.write("<!-- end of EAHPhi_amp -->\n")
        #
        xml_file.write(
            "<!-- EAHPhi_phase: Phase of the effective antenna height in e_phi direction -->\n"
        )
        for f in freqList:
            xml_file.write('<EAHPhi_phase unit="deg" idfreq="%.1f">' % f)
            xml_file.write(nparray2string(vel_ephi_phase[f]) + "</EAHPhi_phase>\n")
        xml_file.write("<!-- end of EAHPhi_phase -->\n")

    print(
        '[INFO] XML file with Isotropic antenna pattern at "'
        + antennaPatternXMLfile
        + '" should have been created.'
    )
    # You should have the XML with the antenna pattern created by now! Check the folder.


if __name__ == "__main__":

    ap = argparse.ArgumentParser()
    ap._action_groups.pop()
    optional = ap.add_argument_group("optional arguments")

    optional.add_argument(
        "-o",
        "--outputpath",
        nargs="?",
        default="../antenna_setup_files/Isotropic_antenna_pattern.xml",
        help="Save path",
    )

    args = vars(ap.parse_args())
    antennaPatternXMLfile = args["outputpath"]
    sys.exit(main())
