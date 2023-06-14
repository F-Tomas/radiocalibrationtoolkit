#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tabulate ULSA (Ultra Longitudinal Sky Absorption) sky maps for a given frequency range and index type.

Usage:
  python <script_name>.py -i <index_type> -f <from_freq> <to_freq> [-s <freq_step>] [-n <nside>] [-o <output_path>]

Arguments:
  -i, --indextype        Type of index used for sky map generation. Possible values are:
                           freq_dependent_index: frequency-dependent spectral index
                           constant_index: constant spectral index
                           direction_dependent_index: direction-dependent spectral index
  -f, --frequencyrange   Frequency range in MHz for which sky maps are to be generated. Should be specified as 'from' and 'to' values.
  -s, --frequencystep    Frequency step for generating sky maps. Default is 1 MHz.
  -n, --nside            Value of the NSIDE parameter for the healpy map. Default is 64.
  -o, --outputpath       Output folder path for saving the generated ULSA sky maps. If not specified, the maps will be saved in the current working directory.

Returns:
  pandas DataFrame: DataFrame containing the ULSA sky maps for the given frequency range and index type. The sky maps are tabulated and saved as a CSV file with filename format: ULSA_<index_type>_maps_<from_freq>-<to_freq>MHz.csv.

Dependencies:
  - ULSA.sky_map.produce_absorbed_sky_map: Module for generating ULSA sky maps
  - pandas: Library for data manipulation and analysis
  - numpy: Library for numerical computing
  - tqdm: Library for progress bars

Example:
  To generate ULSA sky maps for a frequency range of 100-150 MHz with constant spectral index and save the maps in the folder '/maps', run the following command:
  python <script_name>.py -i constant_index -f 100 150 -o /maps
"""

from ULSA.sky_map.produce_absorbed_sky_map import absorption_JRZ
import pandas as pd
import numpy as np
from tqdm import tqdm

import sys
import os
import argparse
from pathlib import Path

import decimal

# Disable
def _blockPrint():
    sys.stdout = open(os.devnull, "w")


# Restore
def _enablePrint():
    sys.stdout = sys.__stdout__


def _frange(i, f, s, endpoint=True):
    d = decimal.Decimal(str(s))
    exp = float(d.as_tuple().exponent)
    print(exp)
    i = int(i * 10**-exp)
    f = int(f * 10**-exp)
    s = int(s * 10**-exp)

    if endpoint:
        f = f + s
    values = range(i, f, s)
    if exp == 0:
        return values
    else:
        return [x * 10**exp for x in values]


def tabulate_ulsa(freq_range, index_type, nside=64):

    ultralon_maps_dict = {}

    for f in tqdm(freq_range):
        _blockPrint()  # supress ULSA output
        cla = absorption_JRZ(
            f,
            nside=nside,
            index_type=index_type,
            distance=50,
            using_raw_diffuse=False,
            using_default_params=True,
            critical_dis=False,
            output_absorp_free_skymap=False,
        )
        ultralon_maps_dict[f] = np.asarray(cla.mpi())[:, 1]
    return pd.DataFrame(ultralon_maps_dict)


def main():

    ultralon_maps_DF = tabulate_ulsa(freq_range, index_type, nside=nside)
    ultralon_maps_DF.to_csv(
        os.path.join(output_path, "ulsa_{}_{:.0f}-{:.0f}MHz.csv").format(
            index_type, freq_range[0], freq_range[-1]
        ),
        index=False,
    )


if __name__ == "__main__":

    ap = argparse.ArgumentParser()
    ap._action_groups.pop()
    required = ap.add_argument_group("required arguments")
    optional = ap.add_argument_group("optional arguments")

    required.add_argument(
        "-i",
        "--indextype",
        nargs="?",
        required=True,
        help="freq_dependent_index constant_index direction_dependent_index ",
    )
    required.add_argument(
        "-f",
        "--frequencyrange",
        nargs="*",
        required=True,
        help="Frequency range in MHz, 'from' 'to'",
    )
    optional.add_argument(
        "-s",
        "--frequencystep",
        nargs="?",
        default=1,
        help="Frequency step. Default is 1 MHz",
    )
    optional.add_argument(
        "-n",
        "--nside",
        nargs="?",
        default=64,
        help="Value of the NSIDE parameter for the healpy map. Default is 64.",
    )
    optional.add_argument(
        "-o",
        "--outputpath",
        nargs="?",
        default="./",
        help="Output folder. If it does not exists, it will be created.",
    )
    args = vars(ap.parse_args())

    index_types = [
        "freq_dependent_index",
        "constant_index",
        "direction_dependent_index",
    ]
    index_type = args["indextype"]
    if index_type not in index_types:
        print("[ERROR] Index type not recognized, use one of these: ", index_types)
        sys.exit()

    freq_range = _frange(
        float(args["frequencyrange"][0]),
        float(args["frequencyrange"][1]),
        float(args["frequencystep"]),
    )
    nside = int(args["nside"])
    output_path = args["outputpath"]

    sys.exit(main())
