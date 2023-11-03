#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  convertSParameterFile2XML.py

"""
convertSParameterFile2XML.py

This script converts S-parameter data from a file into an XML format. It calculates complex voltage amplification, magnitude in dB, and phase in degrees or radians and saves the results in an XML file.

Parameters
----------
- sfile (str): Path to the S-parameter file.
- output_name (str, optional): File name of the produced XML file and the name of the HW response. By default, the filename of the S-parameter file will be used.
- output_dir (str, optional): Output directory. Default is the current directory.
- unwrap (bool, optional): If True, phase values are unwrapped; otherwise, they are left wrapped. Default is True.
- method (str, optional): The amplification method, either 'withS11' (S21/(S11+1)) or 'withoutS11' (S21). Default is 'withS11'.
- rounding (int, optional): Rounding of saved numbers. Default is 4. Use None for no rounding.

Returns
-------
int: Return code, 0 on success.

Note
----
This script requires the `pandas`, `numpy`, `matplotlib`, `scipy`, and `scikit-rf` (skrf) libraries.

Example Usage
-------------
$ python convertSParameterFile2XML.py -s input.s2p -o output.xml -d /output/directory -u True -m withS11 -r 3
"""

import argparse
import re
import io
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import skrf as rf
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit



def main(args):
    np.set_printoptions(suppress=True)
    skrf = rf.Network(sparameterFilePath)
    MHz = 1e-6
    frequencyMHz = skrf.f * MHz

    if method == 'withS11':
        complex_voltage_amplification = (skrf.s21.s/(1+skrf.s11.s)).flatten()
    elif method == 'withoutS11':
        complex_voltage_amplification = (skrf.s21.s/(1)).flatten()
    dBMagnitude = 20*np.log10(np.abs(complex_voltage_amplification))
    
    phaseDegreeWrapped = np.rad2deg(np.angle(complex_voltage_amplification))

    # unwrape phase if wraped
    phaseDegreeUnwrapped = np.unwrap(phaseDegreeWrapped, 180)

    if useUnwrappedPhase == True:
        phaseDeg = phaseDegreeUnwrapped
    else:
        phaseDeg = phaseDegreeWrapped

    #if rounding != None:
     #   dBMagnitude = np.asarray([round(i, rounding) for i in dBMagnitude])
   #     phaseDeg = np.asarray([round(i, rounding) for i in phaseDeg])

    # save in the Offline XML format
    freqString = re.sub(
        " +",
        " ",
        np.array2string(
            frequencyMHz, separator=" ", max_line_width=np.inf, suppress_small=False, formatter={'float_kind': lambda x: f'{x:.{rounding}f}'}
        ).strip("[,]"),
    )
    magString = re.sub(
        " +",
        " ",
        np.array2string(
            dBMagnitude, separator=" ", max_line_width=np.inf, suppress_small=False, formatter={'float_kind': lambda x: f'{x:.{rounding}f}'}
        ).strip("[,]"),
    )
    phaseString = re.sub(
        " +",
        " ",
        np.array2string(
            phaseDeg, separator=" ", max_line_width=np.inf, suppress_small=False,  formatter={'float_kind': lambda x: f'{x:.{rounding}f}'}
        ).strip("[,]"),
    )

    # header = "<!-- {} -->".format(responseName)
    header = ''
    
    responseId = responseName

    freqUnit = "MHz"  # 'MHz'
    magUnit = (
        "dB"
    )  # 'LgAmp'  dB::: <!-- amplification in dB -->  LgAmp::: <!-- amplification in Log10(Amplitude) -->
    phaseUnit = "degree"  #'degree' # 'radiands' ?

    string = (
        header
        + "\n\n\n"
        + '  <RResponse id="'
        + responseId
        + '"> \n'
        + "    <Response> \n "
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
        + '      <Phase unit ="'
        + phaseUnit
        + '">\n'
        "        "
        + phaseString
        + "      </Phase>\n"
        + "    </Response>\n"
        + "  </RResponse>"
    )

    with open(os.path.join(outdir, responseName + ".xml"), "w") as file:
        file.write(string)

    # control plots
    '''
    plt.figure(0)
    skrf.plot_s_db()
    plt.figure(1)
    skrf.plot_s_deg()

    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    fig.suptitle("what goes to the converted XML file")
    ax[0].plot(frequencyMHz, dBMagnitude)
    ax[1].plot(frequencyMHz, phaseDeg)
    [a.set_xlabel("frequency [MHz]") for a in ax]
    ax[0].set_ylabel("dB")
    ax[1].set_ylabel("phase [deg]")
    plt.subplots_adjust(wspace=0.4)
    plt.show()
    '''
    return 0


if __name__ == "__main__":
    import sys
    import os

    # Construct the argument parser
    ap = argparse.ArgumentParser()
    ap._action_groups.pop()
    required = ap.add_argument_group("required arguments")
    optional = ap.add_argument_group("optional arguments")

    # Add the arguments to the parser
    required.add_argument(
        "-s", "--sfile", required=True, help="Path to the s-parameter file."
    )
    optional.add_argument(
        "-o",
        "--output_name",
        nargs="?",
        required=False,
        default=None,
        help="File name of the produced XML file. Also it is the name of the HW response. By default filename of the s-parameter file will be used.",
    )
    optional.add_argument(
        "-d",
        "--output_dir",
        nargs="?",
        required=False,
        default='./',
        help="Output directory",
    )
    optional.add_argument(
        "-u",
        "--unwrap",
        nargs="?",
        required=False,
        default=True,
        help="Unwrap phase. Default is True.",
    )
    optional.add_argument(
        "-m",
        "--method",
        nargs="?",
        required=False,
        default='withS11',
        help="The amplification is either given as S21/(S11+1) or as just S21. Thus our options are 'withS11' or 'withoutS11. Default: 'withS11",
    )
    optional.add_argument(
        "-r",
        "--rounding",
        nargs="?",
        required=False,
        default=4,
        help="Rouding of the save numbers. Default is 4. For no rouding type None",
    )
    args = vars(ap.parse_args())

    method = args['method']
    sparameterFilePath = args["sfile"]
    useUnwrappedPhase = args["unwrap"]
    rounding = int(args["rounding"])

    if args["output_name"] == None:
        responseName = os.path.splitext(os.path.basename(sparameterFilePath))[0]
    else:
        responseName = args["output_name"]

    outdir = args["output_dir"]
    # workDir ='/home/tomas/Documents/scripts/pyUtilScripts/'
    # inputFileName = 'AugerRD_LNA_GRF4002_NS_39mA_ReIm.s2p'
    # sparameterFilePath = workDir+inputFileName
    # useUnwrappedPhase = True
    # rounding = 3

    sys.exit(main(sys.argv))
