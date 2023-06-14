#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse

# This uses the relative path to the package, one directory up.
# Add path to the radiocalibrationtoolkit to your PYTHONPATH to avoid this.
sys.path.append(os.path.join(os.path.abspath(""), "../.."))
from radiocalibrationtoolkit import *


def save_figure(df, ylabel):
    # save figure
    fig = px.imshow(
        df.T,
        width=600,
        aspect="cube",
        color_continuous_scale="jet",
        title="<b>{}</b>".format(gal_model),
    )
    fig.update_layout(
        xaxis=dict(title="<b>LST</b>", tickprefix="<b>", ticksuffix="</b>", dtick=2),
        yaxis=dict(
            title="<b>frequency [MHz]</b>",
            tickprefix="<b>",
            ticksuffix="</b>",
            range=(30, 80),
            tick0=0,
            dtick=10,
            autorange=False,
        ),
    )
    fig.update_layout(
        coloraxis=dict(
            colorbar=dict(
                title=dict(
                    text="<b>{}</b>".format(ylabel),
                    side="right",
                ),
                tickprefix="<b>",
                ticksuffix="</b>",
            ),
            cmin=0,
            cmax=20,
        ),
        font=dict(
            # family=font,
            size=16,
            color="black",
        ),
    )
    fig.write_image(
        os.path.join(
            args["outputpath"]
            + "{}_{}_{}.png".format(
                args["datasettype"], args["savefilenameprefix"], gal_model
            )
        )
    )


def main():
    antenna_inst = AntennaPattern(args["antennamodel"])

    freq_Mhz_range = range(args["frequencyrange"][0], args["frequencyrange"][1], 1)

    if args["datasettype"] == "power":
        lst_range = np.asarray(list(range(lst_range_arg))) + shift_LST_bins

        # read hw file
        hw_file_path = args["hwmodel"]
        hw_dict = read_hw_file(hw_file_path)
        impedance_func = hw_dict["IImpedance"]["antenna_{}".format(orientation)]

        power_density_DF = calculate_power_spectral_density(
            antenna_inst=antenna_inst,
            galactic_map_inst=galactic_map_inst,
            lst_range=lst_range,
            freq_Mhz_range=freq_Mhz_range,
            latitude=args["latitude"],
            update_antenna_conventions={
                "shift_phi": shift_phi,
                "flip_theta": True,
                "flip_phi": False,
                "in_degrees": True,
                "add_invisible_sky": True,
            },
            impedance_func=impedance_func,
        )

        power_DF = integrate_spectral_density(
            power_density_DF, integrated_MHz_bands=power_density_DF.columns
        )

        # apply HW response
        hw_reponse_1 = dB2PowerAmp(hw_dict["RResponse"]["LNA"](power_DF.columns))
        hw_reponse_2 = dB2PowerAmp(hw_dict["RResponse"]["digitizer"](power_DF.columns))
        hw_reponse_3 = dB2PowerAmp(
            hw_dict["RResponse"]["cable_fromLNA2digitizer"](power_DF.columns)
        )
        hw_reponse_4 = dB2PowerAmp(
            hw_dict["RResponse"]["impedance_matching_{}".format(orientation)](
                power_DF.columns
            )
        )

        power_DF = power_DF.multiply(
            hw_reponse_1 * hw_reponse_2 * hw_reponse_3 * hw_reponse_4
        )

        # to piko
        power_DF = power_DF * 1e12

        # save path
        save_path = os.path.join(
            args["outputpath"], "{}_{}".format(args["savefilenameprefix"], gal_model)
        )
        (power_DF.round(3)).to_csv(save_path + ".csv")
        save_figure(power_DF, "Power [pW]")

    elif args["datasettype"] == "voltage2":
        lst_range = range(25)

        voltage2_density_DF = calculate_power_spectral_density(
            antenna_inst=antenna_inst,
            galactic_map_inst=galactic_map_inst,
            lst_range=lst_range,
            freq_Mhz_range=freq_Mhz_range,
            latitude=args["latitude"],
            update_antenna_conventions={
                "shift_phi": shift_phi,
                "flip_theta": True,
                "flip_phi": False,
                "in_degrees": True,
                "add_invisible_sky": True,
            },
        )
        voltage2_density_DF.index.name = "LST"
        voltage2_density_DF.columns.name = "freq_MHz"

        voltage2_density_DF.to_csv(
            os.path.join(
                args["outputpath"],
                "voltage2_density_{}_{}.csv".format(
                    args["savefilenameprefix"], gal_model
                ),
            )
        )


if __name__ == "__main__":
    # select galactic radio emission model
    galactic_models = [
        "LFSS",
        "GSM08",
        "GSM16",
        "Haslam",
        "LFmap",
        "SSM",
        "GMOSS",
        "ULSA",
    ]

    ap = argparse.ArgumentParser()
    ap._action_groups.pop()
    required = ap.add_argument_group("required arguments")
    optional = ap.add_argument_group("optional arguments")

    required.add_argument(
        "-g",
        "--galacticmodel",
        nargs="?",
        required=True,
        choices=galactic_models,
        help="Galactic radio emission model, choose from these: {}".format(
            galactic_models
        ),
    )
    required.add_argument(
        "-a",
        "--antennamodel",
        nargs="?",
        required=True,
        help="Path to the antenna model file.",
    ),
    optional.add_argument(
        "-w",
        "--hwmodel",
        nargs="?",
        default="",
        help="Path to the hardware response file.",
    ),
    optional.add_argument(
        "-f",
        "--frequencyrange",
        nargs="+",
        default=[30, 81],
        type=int,
        help="Frequency range in MHz, 'from' 'to', must be an integer! Default is [30 81)",
    )
    optional.add_argument(
        "-l",
        "--latitude",
        nargs="?",
        default=-35.206667,
        type=float,
        help="Latitude of the local observer, default it -35.206667 (Malargue, Argentina).",
    )
    optional.add_argument(
        "-o",
        "--outputpath",
        nargs="?",
        default="./simulated_power_datasets/",
        help="Output folder. If it does not exists, it will be created.",
    )
    optional.add_argument(
        "-s",
        "--savefilenameprefix",
        nargs="?",
        default="power_",
        help="Saved power dataset filename prefix.",
    )
    optional.add_argument(
        "-t",
        "--datasettype",
        nargs="?",
        default="power",
        help="Calculate 'power' or 'voltage2'.",
    )
    optional.add_argument(
        "-r",
        "--orientation",
        nargs="?",
        default="EW",
        help="This will read HW response with suffix either EW or NS. Default: EW",
    )
    optional.add_argument(
        "-z",
        "--shiftazimuth",
        nargs="?",
        default=-90,
        type=float,
        help="Modify antenna response conventions by shifting the azimuth. Default:-90 degrees.",
    )
    optional.add_argument(
        "-ls",
        "--shiftlst",
        nargs="?",
        default=0.5,
        type=float,
        help="Shift LST bins.",
    )
    optional.add_argument(
        "-lr",
        "--lstrange",
        nargs="?",
        default=24,
        type=int,
        help="Usually one need just 0-23 (or 0.5-23.5) LST. But if the dataset is to be interpolated later, better is to do 0-24 (set the arg to 25 for this)",
    )
    args = vars(ap.parse_args())

    if (args["datasettype"] == "voltage2") and (
        args["outputpath"] == "./simulated_power_datasets/"
    ):
        args["outputpath"] = "./voltage2_density/"

    orientation = args["orientation"]
    gal_model = args["galacticmodel"]
    shift_phi = args["shiftazimuth"]
    shift_LST_bins = args["shiftlst"]
    lst_range_arg = args["lstrange"]

    if gal_model == "LFSS":
        galactic_map_inst = LowFrequencySkyModel(freq_unit="MHz")
    elif gal_model == "GSM08":
        galactic_map_inst = GlobalSkyModel(freq_unit="MHz")
    elif gal_model == "GSM16":
        galactic_map_inst = GlobalSkyModel2016(freq_unit="MHz")
    elif gal_model == "Haslam":
        galactic_map_inst = HaslamSkyModel(freq_unit="MHz", spectral_index=-2.53)
    elif gal_model == "LFmap":
        galactic_map_inst = LFmap()
    elif gal_model == "SSM":
        galactic_map_inst = SSM()
    elif gal_model == "GMOSS":
        galactic_map_inst = GMOSS()
    elif gal_model == "ULSA":
        galactic_map_inst = ULSA(index_type="freq_dependent_index")
    else:
        print("[ERROR] Model not recognized, use some from these:", galactic_models)
        sys.exit()

    mkdir(args["outputpath"])

    sys.exit(main())
