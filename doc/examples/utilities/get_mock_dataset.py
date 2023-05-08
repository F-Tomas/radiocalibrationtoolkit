#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from radiocalibrationtoolkit import *
import argparse
from pathlib import Path

piko = 1e-12

def to_string(input_val):
    if isinstance(input_val, (float, int)):
        return str(input_val)
    elif isinstance(input_val, Union[tuple, list]):
        return "_".join(str(x) for x in input_val)
    else:
        raise TypeError("Input must be a float, int, or tuple")


def main():
    # read HW response

    hw_dict = read_hw_file(
        args["hwprofilefile"], interp_args={"fill_value": "extrapolate"}
    )
    print()

    hw_reponse_1 = hw_dict["RResponse"]["LNA"]
    hw_reponse_2 = hw_dict["RResponse"]["digitizer"]
    hw_reponse_3 = hw_dict["RResponse"]["cable_fromLNA2digitizer"]
    hw_reponse_4 = hw_dict["RResponse"]["impedance_matching_EW"]

    # merge all hw responses to one function
    def hw_response_func(x):
        return dB2PowerAmp(
            hw_reponse_1(x) + hw_reponse_2(x) + hw_reponse_3(x) + hw_reponse_4(x)
        )

    # impedance function
    impedance_func = hw_dict["IImpedance"]["antenna_EW"]

    # read sidereal voltage square spectral density
    sidereal_voltage2_density_DF = pd.read_csv(
        args["voltagedensity2file"],
        index_col=0,
    )
    sidereal_voltage2_density_DF.columns = sidereal_voltage2_density_DF.columns.astype(
        float
    )

    def hw_response_func(x):
        return dB2PowerAmp(
            hw_reponse_1(x) + hw_reponse_2(x) + hw_reponse_3(x) + hw_reponse_4(x)
        )

    mock_trace_generator = Mock_trace_generator(
        sidereal_voltage2_density_DF=sidereal_voltage2_density_DF,
        hw_response_func=hw_response_func,
        impedance_func=impedance_func,
        voltage2ADC=2048,
        time_trace_size=2048,
        sampling_frequency_MHz=250,
    )
    freq_MHz_bins = mock_trace_generator.get_frequency_bins()

    mock_traces_DF = mock_trace_generator.generate_mock_trace(
        args["numberoftraces"],
        temp_celsius=args["temperature"],
        additional_noise=args["additionalnoise"],
        # nbi={"67.25": 1},
        # nbi_err=0.3,
        rounding=not args["turn_off_rounding"],
    )

    # system parameters
    sampling_frequency_MHz = 250
    N = mock_traces_DF.columns[2:].size
    ADC2Volts = 1 / 2048
    trace_time_length_sec = N / (sampling_frequency_MHz * 1e6)

    # FFT to spectra
    spectra_df = time_trace_df_2_spectra_df(
        mock_traces_DF, DSC=2, sampling_frequency_MHz=sampling_frequency_MHz
    )

    # use the formula, create the integrand first
    integrand_df = (
        (spectra_df * ADC2Volts / (sampling_frequency_MHz * 1e6)) ** 2
    ).divide(impedance_func(spectra_df.columns.values))

    fr = args["frequencyrange"]
    # integrate
    mock_power_unbinned_DF = (2 / trace_time_length_sec) * integrate_spectral_density(
        integrand_df,
        integrated_MHz_bands=np.linspace(fr[0], fr[1], fr[1] - fr[0] + 1),
        integrating_method="on_discontinuous_function",
    )

    mock_power_unbinned_DF = pd.concat(
        (mock_traces_DF.iloc[:, :2], mock_power_unbinned_DF), axis=1
    )

    mock_power_DF = bin_df_rows(
        mock_power_unbinned_DF, binning_column="lst", bins=list(range(25))
    )
    mock_power_DF.index.name = "lst"
    mock_power_DF = mock_power_DF.drop(["temp_c", "lst"], axis=1)

    if args["overridesavefilename"] is not None:
        filename = args["savefilename"]
    else:
        filename = (
            "mock_power_dataset-{}_N{}_temp{}C_{}additionalnoise_rounding-{}.csv".format(
                Path(args["voltagedensity2file"]).stem.replace("voltage2_density_", ""),
                args["numberoftraces"],
                to_string(args["temperature"]),
                args["additionalnoise"],
                not args["turn_off_rounding"],
            )
        )

    ((mock_power_DF/piko).round(3)).to_csv(os.path.join(args["outputpath"], filename))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap._action_groups.pop()
    required = ap.add_argument_group("required arguments")
    optional = ap.add_argument_group("optional arguments")

    required.add_argument(
        "-v",
        "--voltagedensity2file",
        nargs="?",
        required=True,
        help="Path to the antenna model file.",
    ),
    optional.add_argument(
        "-w",
        "--hwprofilefile",
        nargs="?",
        default="",
        help="Path to the hardware response file.",
    ),
    optional.add_argument(
        "-f",
        "--frequencyrange",
        nargs="*",
        default=(30, 81),
        help="Frequency range in MHz, 'from' 'to', must be an integer! Default is [30 81)",
    )
    optional.add_argument(
        "-t",
        "--temperature",
        nargs="+",
        default=[-10, 50],
        action=AppendOnlyIfMore,
        type=float,
        help="Temperature range. Random temperature from this range will be used for the temperature noise.\
                If only one values is specified,\
                the temperature will be fixed to this temperature.",
    )
    optional.add_argument(
        "-n",
        "--numberoftraces",
        nargs="?",
        default=5000,
        type=int,
        help="Number of generated time traces",
    )
    optional.add_argument(
        "-a",
        "--additionalnoise",
        nargs="?",
        default=0,
        type=float,
        help="Mean of the additional noise added to final spectrum. Default is zero (no extra noise).\
            Reasonable value for the extra noise is 5e-16",
    )
    optional.add_argument(
        "-tor",
        "--turn-off-rounding",
        # nargs="?",
        action="store_true",
        help="Whether the generated time traces in ADC units should be rounded.",
    )
    optional.add_argument(
        "-o",
        "--outputpath",
        nargs="?",
        default="./mock_power_datasets/",
        help="Output folder. If it does not exists, it will be created.",
    )
    optional.add_argument(
        "-s",
        "--overridesavefilename",
        nargs="?",
        default=None,
        help="Overrides filename.",
    )

    args = vars(ap.parse_args())
    mkdir(args["outputpath"])

    print(args)
    sys.exit(main())
