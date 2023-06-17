#!/usr/bin/env python
# -*- coding: utf-8 -*-

## import sys
## import os
## _path = os.path.dirname(os.path.abspath(__file__))
## sys.path.append(os.path.join(_path, '../../..')
## from radiocalibrationtoolkit import *

from scipy import constants
from scipy import integrate
import statsmodels.api as sm
from typing import Union
from tqdm import tqdm
from ..antennareader import *
from ..skymapwrappers import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator

NSIDE = 64
pixel_theta, pixel_phi = hp.pix2ang(NSIDE, np.arange(hp.nside2npix(NSIDE)))
ANGULAR_FACTOR = hp.pixelfunc.nside2pixarea(NSIDE) * np.sin(pixel_theta)
CONSTANT_FACTOR = (2 * constants.k / constants.c**2) * constants.physical_constants[
    "characteristic impedance of vacuum"
][0]


def voltage_squared_spectral_density(
    antenna_map: np.array, galactic_map: np.array, frequency_MHz: float
) -> float:
    """
    Calculate the voltage squared spectral density for a given antenna map and galactic map at a specific frequency.

    Parameters
    ----------
    antenna_map : np.array
        Array representing the antenna map.
    galactic_map : np.array
        Array representing the galactic map.
    frequency_MHz : float
        Frequency in MHz.

    Returns
    -------
    float
        Result of the voltage squared spectral density calculation.

    """
    return (
        (frequency_MHz * 1e6) ** 2
        * CONSTANT_FACTOR
        * integrate_hpmap(galactic_map * antenna_map**2)
    )


def add_integral_edges_to_x_axis(
    x: Union[list, np.ndarray], integral_edges: Union[list, np.ndarray]
) -> np.ndarray:
    """
    Add integral edges to the x-axis values and return the sorted unique array.

    Parameters
    ----------
    x : list or np.ndarray
        Input array representing the x-axis values.
    integral_edges : list or np.ndarray
        Integral edges to be added to the x-axis values.

    Returns
    -------
    np.ndarray
        Sorted unique array of the x-axis values with integral edges.

    """
    x_new = np.hstack((x, integral_edges))
    x_new.sort()
    x_new = np.asarray(list(dict.fromkeys(x_new)))
    return x_new


def integrate_on_a_part_of_an_array(
    x: Union[list, np.ndarray],
    interp_func: interp1d,
    lower_integral_boundary: float,
    upper_integral_boundary: float,
) -> float:
    """
    Perform numerical integration on a specific part of an array using the given interpolation function.

    Parameters
    ----------
    x : list or np.ndarray
        Array representing the x-axis values.
    interp_func : interp1d
        Interpolation function that maps x-values to y-values.
    lower_integral_boundary : float
        Lower boundary of the integral.
    upper_integral_boundary : float
        Upper boundary of the integral.

    Returns
    -------
    float
        Result of the numerical integration.

    """
    value_filter = (x >= lower_integral_boundary) & (x <= upper_integral_boundary)
    return integrate.trapz(interp_func(x[value_filter]), x[value_filter])


def integrate_spectral_density(
    sd_DF, integrated_MHz_bands, integrating_method="quad", point=None
):
    """
    Compute the integrated spectral density for each row of a given pandas DataFrame across specified frequency bands.

    Parameters
    ----------
    sd_DF : pandas DataFrame
        DataFrame with spectral density values as the data and frequency bins as the columns.
    integrated_MHz_bands : list or numpy.ndarray
        List or array of floats representing the edges of the frequency bands to integrate over, in units of MHz.
    integrating_method : str, optional
        Integration method to use. Options are "quad" (default) or "on_discontinuous_function".
    point : list, numpy.ndarray or None, optional
        Points to use in the integration method if "quad" is chosen. Defaults to None.

    Returns
    -------
    pandas DataFrame
        DataFrame with the integrated spectral density for each row of the input DataFrame across each frequency band
        specified by the `integrated_MHz_bands` parameter.

    Raises
    ------
    ValueError
        If the input DataFrame has fewer than 2 frequency bins.

    Examples
    --------
    >>> df = pd.DataFrame({'f1': [1, 2, 3], 'f2': [4, 5, 6]})
    >>> integrate_spectral_density(df, [1.5, 4.5])
           1.5       4.5
    0  2.621005  9.453482
    """
    if sd_DF.shape[1] <= 1:
        print(
            "[ERROR] Not enough frequencies to integrate, calculate at least two frequency bins."
        )
    mega = 1e6
    freq_Hz = sd_DF.columns.values * mega
    freq_new_Hz = add_integral_edges_to_x_axis(freq_Hz, integrated_MHz_bands * mega)

    temp = []
    for row in sd_DF.index.values:
        density_func = interp1d(freq_Hz, sd_DF.loc[row, :].values)
        temp_nested = []
        for i, _ in enumerate(integrated_MHz_bands[:-1]):
            if integrating_method == "on_discontinuous_function":
                temp_nested.append(
                    integrate_on_a_part_of_an_array(
                        freq_new_Hz,
                        density_func,
                        integrated_MHz_bands[i] * mega,
                        integrated_MHz_bands[i + 1] * mega,
                    )
                )
            elif integrating_method == "quad":
                # this method is super slow on discontinuous function with 'points'
                temp_nested = temp_nested + [
                    integrate.quad(
                        density_func,
                        integrated_MHz_bands[i] * mega,
                        integrated_MHz_bands[i + 1] * mega,
                        points=point,
                    )[0]
                ]
        temp.append(np.asarray(temp_nested))
    integrated_DF = pd.DataFrame(
        np.asarray(temp),
        index=sd_DF.index,
        columns=integrated_MHz_bands[:-1],
        dtype="float64",
    )
    integrated_DF.index.name = sd_DF.index.name
    integrated_DF.columns.name = "freq_MHz_band"
    return integrated_DF


def calculate_power_spectral_density(
    antenna_inst,
    galactic_map_inst,
    lst_range: list,
    freq_Mhz_range: list,
    latitude: float,
    impedance_func: Union[interp1d, None] = None,
    update_antenna_conventions: dict = {},
    transform_hpmaps2_local_cs: bool = False,
    nside_transform: int = 512,
    distort_antenna_map: float = 0,
    distort_galactic_map: float = 0,
) -> pd.DataFrame:
    """Calculate the power spectral density of the sky signal as received by an antenna.

    Parameters
    ----------
    antenna_inst : Antenna
        Instance of the antenna model used to receive the signal.
    galactic_map_inst : GalacticMap
        Instance of the galactic map model used to generate the sky signal.
    lst_range : List[float]
        List of local sidereal times (in hours) at which to calculate the power spectral density.
    freq_Mhz_range : List[float]
        List of frequencies (in MHz) at which to calculate the power spectral density.
    latitude : float
        Latitude of the observer in degrees.
    impedance_func : Optional[interp1d], default None
        Interpolation function representing the antenna's impedance as a function of frequency.
        If not provided, a unity impedance function is assumed.
    update_antenna_conventions : Dict[str, Union[str, float]], default {}
        Dictionary of update arguments to pass to the `Data2hpmap` constructor when creating the
        antenna map instance. See the `Data2hpmap` documentation for details.
    transform_hpmaps2_local_cs : bool, default False
        If True, transform the galactic map and antenna map to the local coordinate system
        defined by the observer's latitude and the local sidereal time. This approach is not
        recommended due to the high computational cost of transforming the galactic map.
        The default is to use the galactic coordinate system and create a static antenna map.
    nside_transform : int, default 512
        Temporary NSIDE upgrade of the galactic map when using `transform_hpmaps2_local_cs=True`.
        A lower value will increase speed, but the galactic map in local coordinates will have
        higher numerical errors.

    Returns
    -------
    pd.DataFrame
        Dataframe with the calculated power spectral density values. The columns correspond to
        frequencies, and the rows correspond to local sidereal times.
    """
    if impedance_func is None:
        print("[INFO] No impedance function, assuming it is equal to '1'")
        impedance_func = interp1d(freq_Mhz_range, np.ones(len(freq_Mhz_range)))

    power_density_DF = pd.DataFrame(
        index=list(lst_range), columns=freq_Mhz_range, dtype="float64"
    )

    temp = []
    for frequency_MHz in tqdm(freq_Mhz_range):

        galactic_map = hp.ma(
            hp.pixelfunc.ud_grade(galactic_map_inst.generate(frequency_MHz), NSIDE)
        ).copy()
        # these to two aproaches are equivalent, but the second is ~25% faster
        # antenna_map_inst = antenna_inst.convert2hp(frequency=frequency_MHz, quantity='absolute', **update_antenna_conventions)
        antenna_map_inst = Data2hpmap(
            antenna_inst.get(frequency=frequency_MHz, quantity="absolute"),
            **update_antenna_conventions,
        )
        if transform_hpmaps2_local_cs:
            # create the antenna hp map array in local coordinates (not recommened, all maps should be by default in galactic coordinates)
            # now, it this map is NOT a function of LST!
            antenna_map = antenna_map_inst.get_map()
            galactic_map_in_galactic_cs = galactic_map.copy()
        temp_row = []
        for lst in lst_range:
            rotation_parameters = create_rotation_parameters(lst, latitude)
            rotator = create_rotator(lst, latitude, coord=["G", "C"])
            if transform_hpmaps2_local_cs:
                # permanently alter galactic map array to local coordinates
                # now, it this map is a function of LST!
                galactic_map = (
                    rotate_default_hpmap_array_from_galactic2local_coordinates(
                        galactic_map_in_galactic_cs,
                        lst,
                        latitude,
                        nside_transform=nside_transform,
                    )
                )
            else:
                galactic_map.mask = create_local_mask(NSIDE, rotation_parameters)
                antenna_map = antenna_map_inst.get_map(rotator=rotator)
            # todo fragmentation warn
            #  power_density_DF.loc[lst, frequency_MHz] =
            temp_row.append(
                voltage_squared_spectral_density(
                    antenna_map=distort_array(antenna_map, rstd=distort_antenna_map),
                    galactic_map=distort_array(galactic_map, rstd=distort_galactic_map),
                    frequency_MHz=frequency_MHz,
                )
                / impedance_func(frequency_MHz)
            )
        temp.append(np.asarray(temp_row))
    power_density_DF.iloc[:, :] = np.asarray(temp).T
    return power_density_DF


def time_trace_df_2_spectra_df(time_trace_df, DSC=0, sampling_frequency_MHz=250):
    # DSC - data start column, at which column the time trace starts
    """
    Converts a time trace DataFrame to a spectral density DataFrame using FFT.

    Parameters:
    -----------
    time_trace_df : pandas.DataFrame
        DataFrame containing the time trace data. The signal should be stored in columns.
    DSC : int, optional (default = 0)
        The data start column, i.e., the index of the first column containing the signal.
    sampling_frequency_MHz : float, optional (default = 250)
        Sampling frequency of the signal in MHz.

    Returns:
    --------
    pandas.DataFrame
        DataFrame containing the spectral density of the signal. The frequency values are stored in the columns.
    """
    N = time_trace_df.columns[DSC:].size
    spectra_df = pd.DataFrame(
        [i for i in time_trace_df.iloc[:, DSC:].apply(abs_fft, axis=1)]
    )
    spectra_df.columns = np.arange(N // 2 + 1) / (N / 2) * sampling_frequency_MHz / 2
    return spectra_df


def get_bootstrap_CI(data, B=5000, alpha=1 - 0.6826, statistic_func=np.median):
    """
    Calculates the bootstrap confidence interval for the median of the given data.

    Parameters:
    -----------
    data : numpy.ndarray
        One-dimensional array containing the data.
    B : int, optional (default = 1000)
        The number of bootstrap samples.
    alpha : float, optional (default = 0.32)
        The level of significance, i.e., the proportion of the confidence interval to be calculated.

    Returns:
    --------
    tuple
        Lower and upper bounds of the confidence interval.
    """
    median_values = np.zeros(B)
    for i in range(B):
        sample = np.random.choice(data, size=len(data), replace=True)
        central_values[i] = statistic_func(sample)

    # Calculate confidence interval
    lower_percentile = alpha / 2.0
    upper_percentile = 1 - lower_percentile
    lower_bound = np.percentile(central_values, 100 * lower_percentile)
    upper_bound = np.percentile(central_values, 100 * upper_percentile)
    return lower_bound, upper_bound


def robust_regression(x_arr, y_arr):
    """
    Calculates the robust regression coefficients for the given data using the HuberT norm.

    Parameters:
    -----------
    x_arr : numpy.ndarray
        One-dimensional array containing the predictor data.
    y_arr : numpy.ndarray
        One-dimensional array containing the response data.

    Returns:
    --------
    numpy.ndarray
        Array containing the intercept and slope of the regression model.
    """
    rlm = sm.RLM(y_arr, sm.add_constant(x_arr), M=sm.robust.norms.HuberT())
    rlm_results = rlm.fit()
    return rlm_results.params  # intercept, slope


def dB2PowerAmp(dB):
    """
    Converts the given value from decibel scale to power amplitude scale.

    Parameters:
    -----------
    dB : float
        Value to be converted in decibel scale.

    Returns:
    --------
    float
        Value converted to power amplitude scale.
    """
    return 10 ** (dB / 10)


def dB2VoltageAmp(dB):
    """
    Converts the given value from decibel scale to voltage amplitude scale.

    Parameters:
    -----------
    dB : float
        Value to be converted in decibel scale.

    Returns:
    --------
    float
        Value converted to voltage amplitude scale.
    """
    return 10 ** (dB / 20)


def powerAmp2dB(powerAmp):
    """
    Converts the given value from power amplitude scale to decibel scale.

    Parameters:
    -----------
    powerAmp : float
        Value to be converted in power amplitude scale.

    Returns:
    --------
    float
        Value converted to decibel scale.
    """
    return 10 * np.log10(powerAmp)


def voltageAmp2dB(voltageAmp):
    """
    Converts the given value from voltage amplitude scale to decibel scale.

    Parameters:
    -----------
    voltageAmp : float
        Value to be converted in voltage amplitude scale.

    Returns:
    --------
    float
        Value converted to decibel scale.
    """
    return 20 * np.log10(voltageAmp)


def get_energy_from_two_sided_spectrum(spectrum: np.ndarray) -> float:
    """
    Calculates the energy from a two-sided spectrum.

    Parameters:
        spectrum (np.ndarray): Two-sided spectrum array.

    Returns:
        float: Energy calculated from the spectrum.
    """
    return np.sum(spectrum**2) / spectrum.size


def get_energy_from_time_trace(time_trace: np.ndarray) -> float:
    """
    Calculates the energy from a time trace.

    Parameters:
        time_trace (np.ndarray): Time trace array.

    Returns:
        float: Energy calculated from the time trace.
    """
    return np.sum(time_trace**2)


def get_energy_from_one_sided_spectrum(rspectrum: np.ndarray) -> float:
    """
    Calculates the energy from a one-sided spectrum.

    Parameters:
        rspectrum (np.ndarray): One-sided spectrum array.

    Returns:
        float: Energy calculated from the one-sided spectrum.
    """
    return (
        2 * np.sum(rspectrum[1:-1] ** 2) + rspectrum[0] ** 2 + rspectrum[-1] ** 2
    ) / (2 * (rspectrum.size - 1))


def get_energy_from_one_sided_spectrum_corrected4one_side(
    r2spectrum: np.ndarray,
) -> float:
    """
    Calculates the energy from a one-sided spectrum with correction for one side.

    Parameters:
        r2spectrum (np.ndarray): One-sided spectrum array.

    Returns:
        float: Energy calculated from the one-sided spectrum with correction.
    """
    return np.sum((r2spectrum) ** 2) / (2 * (r2spectrum.size - 1))


def correct_energy_of_one_sided_spectrum(spectrum: np.ndarray) -> np.ndarray:
    """
    Corrects the energy of a one-sided spectrum.

    Parameters:
        spectrum (np.ndarray): One-sided spectrum array.

    Returns:
        np.ndarray: Corrected one-sided spectrum array.
    """
    r2spectrum = spectrum.copy()
    r2spectrum[1:-1] = r2spectrum[1:-1] * np.sqrt(2)
    return r2spectrum


def linspace_with_middle_value(
    middle: float, radius: float, num_samples: int
) -> np.ndarray:
    """
    Generates a 1-D array of evenly spaced values centered around a middle value.

    Parameters:
        middle (float): The middle value.
        radius (float): The radius around the middle value.
        num_samples (int): The number of samples to generate.

    Returns:
        np.ndarray: An array of evenly spaced values.
    """
    return np.linspace(middle - radius, middle + radius, num_samples)


def calculate_truncated_stats(
    data: np.ndarray, lower_percentile: float, upper_percentile: float
) -> tuple:
    """
    Calculates the truncated mean and standard deviation of the given data within the specified percentiles.

    Parameters:
        data (np.ndarray): Input data array.
        lower_percentile (float): Lower percentile threshold.
        upper_percentile (float): Upper percentile threshold.

    Returns:
        tuple: A tuple containing the truncated mean and standard deviation.
    """
    # Filter the data within the specified percentiles
    truncated_data = truncate_data(*locals().values())

    # Calculate the truncated mean and standard deviation
    truncated_mean = np.mean(truncated_data)
    truncated_std = np.std(truncated_data)

    return truncated_mean, truncated_std


def truncate_data(
    data: np.ndarray, lower_percentile: float, upper_percentile: float
) -> np.ndarray:

    # Calculate the lower and upper thresholds
    lower_threshold = np.percentile(data, lower_percentile)
    upper_threshold = np.percentile(data, upper_percentile)

    # Filter the data within the specified percentiles
    return data[(data >= lower_threshold) & (data <= upper_threshold)]


def get_fitted_voltage_calibration_params_and_noise_offsets(
    power_sim_DF: pd.DataFrame, power_rec_DF: pd.DataFrame
) -> tuple:
    """
    Calculate the fitted voltage calibration parameters and noise offsets.

    Parameters:
        power_sim_DF : pd.DataFrame
            Simulated power DataFrame with frequency columns.
        power_rec_DF : pd.DataFrame
            Recorded power DataFrame with frequency columns.

    Returns:
        tuple
            A tuple containing two DataFrames:
            - First DataFrame: Fitted voltage calibration parameters with frequency columns.
            - Second DataFrame: Noise offsets with frequency columns.

    """
    slopes = []
    intercepts = []
    freqs = power_sim_DF.columns
    for i, freq in enumerate(freqs):
        x_arr = power_sim_DF.loc[:, freq].values
        y_arr = power_rec_DF.loc[:, freq].values

        intercept, slope = robust_regression(x_arr, y_arr)
        intercepts.append(intercept)
        slopes.append(slope)

    slopes = np.asarray(slopes) ** (1 / 2)
    intercepts = np.asarray(intercepts)
    return pd.DataFrame([slopes], columns=freqs), pd.DataFrame(
        [intercepts], columns=freqs
    )


def get_fitted_voltage_cal_params_and_noise_offsets_from_concat_sim_dfs(
    concatenated_df: pd.DataFrame, power_rec_DF: pd.DataFrame
) -> tuple:
    """
    Calculate the fitted voltage calibration parameters and noise offsets
    from concatenated simulation and recorded power DataFrames.

    Parameters
    ----------
    concatenated_df : pd.DataFrame
        Concatenated simulation DataFrame with multi-level index.
    power_rec_DF : pd.DataFrame
        Recorded power DataFrame with frequency columns.

    Returns
    -------
    tuple
        A tuple containing two DataFrames:
        - First DataFrame: Fitted voltage calibration parameters with frequency columns.
        - Second DataFrame: Noise offsets with frequency columns.
    """
    slopes_dict = {}
    intercepts_dict = {}

    for key in concatenated_df.index.levels[0]:
        power_sim_DF = concatenated_df.xs(key)
        s, i = get_fitted_voltage_calibration_params_and_noise_offsets(
            power_sim_DF, power_rec_DF
        )
        slopes_dict[key] = s.values.flatten()
        intercepts_dict[key] = i.values.flatten()

    slopes_DF = pd.DataFrame(slopes_dict).T
    slopes_DF.columns = power_sim_DF.columns
    intercepts_DF = pd.DataFrame(intercepts_dict).T
    intercepts_DF.columns = power_sim_DF.columns

    return slopes_DF, intercepts_DF


def get_frequency_independent_calibration_param(
    slopes_DF: pd.DataFrame, show_plots: bool = False
) -> tuple:
    """
    Calculate the frequency-independent calibration parameters.

    Parameters
    ----------
    slopes_DF : pd.DataFrame
        DataFrame containing slopes for each frequency.
    show_plots : bool, optional
        Flag indicating whether to show plots, by default False.

    Returns
    -------
    tuple
        A tuple containing the statistics of the frequency-independent calibration parameters:
        - Total central value
        - Lower bound of the central value
        - Upper bound of the central value

    """
    bounds = (0.2, 2)
    bins = np.linspace(*bounds, 3000)
    trunc = 5
    tolerance = 0.001
    norm_test = 0
    i = 1
    norm_test_list = []
    stats_list = []
    if show_plots:
        fig, ax = plt.subplots()
        ax.set_xlabel("voltage calibration constant")
    while (np.abs(norm_test - 1) > tolerance) & (trunc < 50):
        print("Truncating data to {}, {} percentils".format(trunc, 100 - trunc))
        print("Iteration number: {} out of max 9 iterations".format(i))
        xax_, log_dens, stats = apply_KDE(
            truncate_data(dropnans(slopes_DF.values.flatten()), 5, 100 - trunc),
            bounds=bounds,
            show_plots=show_plots,
        )
        norm_test, _ = quad(
            interp1d(xax_, np.exp(log_dens)),
            bins[0],
            bins[-1],
            epsabs=1e-5,
            epsrel=1e-5,
        )
        norm_test_list.append(norm_test)
        stats_list.append(stats)
        # print(norm_test_list)
        if i > 1:
            if norm_test_list[-1] < norm_test_list[-2]:
                print(
                    "Norm has grown! Breaking the iteration and outputing the statistics from previous iteration!"
                )
                stats = stats_list[-2]
                norm_test = 1  # this is to break the loop
        i += 1
        trunc += 5
    return stats


def get_stats_of_freq_dependent_cal_parameters_from_multiple_sim_datasets(
    slopes_DF: pd.DataFrame,
) -> tuple:
    """
    Calculate the statistics of frequency-dependent calibration parameters from multiple simulation datasets.

    Parameters
    ----------
    slopes_DF : array-like
        The input DataFrame containing slopes for each frequency.

    Returns
    -------
    freq_dependent_slopes : ndarray
        Array of mean slopes for each frequency.
    bounds : ndarray
        Array of standard deviations of slopes for each frequency.

    """
    bounds = []
    freq_dependent_slopes = []
    frequencies_MHz = []
    for freq in slopes_DF.columns:
        bounds.append(np.std(slopes_DF.loc[:, freq].values))
        freq_dependent_slopes.append(np.mean(slopes_DF.loc[:, freq].values))
        frequencies_MHz.append(freq)
    bounds = np.asarray(bounds)
    freq_dependent_slopes = np.asarray(freq_dependent_slopes)
    return freq_dependent_slopes, bounds


def get_and_plot_calibration_results(
    slopes_DF: pd.DataFrame,
    intercepts_DF: pd.DataFrame,
    title: str = "",
    labels: str or list = "",
) -> tuple[float, float, float]:
    """
    Calculate and plot calibration results.

    Parameters
    ----------
    slopes_DF : DataFrame
        DataFrame containing slopes for each frequency.
    intercepts_DF : DataFrame
        DataFrame containing intercepts for each frequency.
    title : str, optional
        Title for the plot (default is an empty string).
    labels : str or list, optional
        Labels for the plot legend. If empty string, labels will be extracted from the index of `slopes_DF`.
        If `None`, no labels will be displayed. (default is an empty string).

    Returns
    -------
    total_central_value : float
        Total central value of the calibration results.
    bounds_total_central_value : float
        Bounds of the total central value.

    """
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    alpha = 0.1
    (
        freq_dependent_slopes,
        bounds,
    ) = get_stats_of_freq_dependent_cal_parameters_from_multiple_sim_datasets(slopes_DF)
    frequencies_MHz = slopes_DF.columns.values.astype(float)

    stats = get_frequency_independent_calibration_param(slopes_DF, show_plots=True)

    total_central_value = stats[2]
    bounds_total_central_value = stats[:2]

    # plot
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(title)

    # for a single map, this does not make sense to do
    if slopes_DF.index.size > 1:
        ax[0].errorbar(
            frequencies_MHz,
            freq_dependent_slopes,
            yerr=bounds,
            label="$\mu(f)$",
            color="k",
        )

    if labels == "":
        glabels = [k.split("_")[2] for k in slopes_DF.index.values]
    elif labels == None:
        glabels = [None] * slopes_DF.index.size
    else:
        glabels = labels

    for i, row in enumerate(slopes_DF.index):
        color = color_cycle[i % len(color_cycle)]
        ax[0].plot(
            frequencies_MHz,
            slopes_DF.loc[row, :].values,
            linestyle="-",
            alpha=alpha,
            color=color,
        )  # , label=glabels)
        ax[0].plot(
            frequencies_MHz,
            slopes_DF.loc[row, :].values,
            marker="o",
            linestyle="",
            markersize=3,
            color=color,
        )  # , label=glabels)

        if slopes_DF.index.size > 1:
            ax[1].plot(
                frequencies_MHz,
                intercepts_DF.loc[row, :].values,
                color=color,
                label=glabels[i],
            )
        else:
            ax[1].bar(
                frequencies_MHz,
                intercepts_DF.loc[row, :].values,
                color=color,
                label=glabels[i],
                width=1,
            )
    ax[0].set_xlabel("frequency [MHz]")
    ax[0].set_ylabel("voltage \ncalibration parameter")
    ax[1].set_ylim(1, 1000)
    ax[1].set_xlabel("frequency [MHz]")
    ax[1].set_ylabel("noise offset [pW]")
    ax[1].set_yscale('log')
    ax[1].xaxis.set_major_locator(MultipleLocator(10))
    if labels is not None:
        ax[1].legend(fontsize=12, ncol=2)

    ax[0].axes.axhline(
        total_central_value,
        color="grey",
        lw=3,
        label=r"$\mu_{{\mathrm{{trun}}}}^{{\mathrm{{KDE}}}}={:.2f}_{{-{:.1f}\%}}^{{+{:.1f}\%}}$".format(
            total_central_value,
            abs(bounds_total_central_value[0] / total_central_value - 1) * 100,
            abs(bounds_total_central_value[1] / total_central_value - 1) * 100,
        ),
    )
    ax[0].axes.axhspan(
        bounds_total_central_value[0],
        bounds_total_central_value[1],
        color="grey",
        alpha=0.3,
        # label="68% CI",
    )
    ax[0].legend(fontsize=12)

    fig.subplots_adjust(
        left=0.2,
        bottom=0.2,
        wspace=0.5,
    )

    return total_central_value, *bounds_total_central_value


def apply_KDE(
    input_data: pd.DataFrame,
    bounds: list | tuple | np.ndarray = [0.5, 1.5],
    show_plots: bool = True,
    data_in_relative_values: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply Kernel Density Estimation (KDE) on input data and calculate relevant statistics.

    Parameters:
    -----------
    input_data : pd.DataFrame
        Input data for KDE estimation.

    bounds : np.ndarray, optional
        Bounds of the KDE estimation range. Default is (0.5, 1.5).

    show_plots : bool, optional
        Flag to indicate whether to show the resulting plots. Default is True.

    data_in_relative_values : bool, optional
        Flag to indicate whether the input data is in relative values. Default is False.

    Returns:
    --------
    tuple[np.ndarray, np.ndarray]
        Tuple containing the x-axis values, log densities, and calculated statistics.

    Notes:
    ------
    - This function applies Kernel Density Estimation (KDE) on the input data using a Gaussian kernel and bandwidth of 0.05.
    - It calculates the log densities and relevant statistics based on the KDE estimation.
    - If show_plots is True, it displays the KDE plot along with relevant annotations and fills the area within +/- one standard deviation.

    Example:
    --------
    # Example usage
    xaxis_values, log_densities, statistics = apply_KDE(input_data, bounds=(0.2, 1.8), show_plots=True)
    """
    xax_kde = np.linspace(bounds[0], bounds[-1], 2000)[:, np.newaxis]
    kde = KernelDensity(kernel="gaussian", bandwidth=0.05).fit(
        input_data[:, np.newaxis]
    )
    log_dens = kde.score_samples(xax_kde)

    xax_ = xax_kde.flatten()

    center_i = find_index_on_CDF(log_dens, xax_kde, 0.5)
    left_i = find_index_on_CDF(log_dens, xax_kde, 0.5 - 0.341)
    right_i = find_index_on_CDF(log_dens, xax_kde, 0.5 + 0.341)

    if data_in_relative_values:
        label = r"$\mu={:.2f}_{{-{:.2f}}}^{{+{:.2f}}}$".format(
            xax_[center_i],
            abs(xax_[left_i + 1] - xax_[center_i]),
            abs(xax_[right_i + 1] - xax_[center_i]),
        )
    else:
        label = r"$\mu={:.2f}_{{-{:.2f}\%}}^{{+{:.2f}\%}}$".format(
            xax_[center_i],
            abs(xax_[left_i + 1] / xax_[center_i] - 1) * 100,
            abs(xax_[right_i + 1] / xax_[center_i] - 1) * 100,
        )
    stats = (xax_[left_i + 1], xax_[right_i + 1], xax_[center_i])

    print(
        "Test1: p0.5={}".format(
            np.sum(np.exp(log_dens[: center_i + 1]) * np.diff(xax_)[0]) - 0.5
        )
    )
    print(
        "Test2: p0.5={}".format(
            np.sum(np.exp(log_dens[center_i + 1 :]) * np.diff(xax_)[0]) - 0.5
        )
    )
    print(
        "Test3: p0.5+/-0.341={}".format(
            np.sum(np.exp(log_dens[left_i + 1 : right_i + 1]) * np.diff(xax_)[0])
            - 2 * 0.341
        )
    )
    print(xax_[center_i], xax_[left_i + 1], xax_[right_i])

    if show_plots:
        fig = plt.gcf()
        if not fig.axes:
            ax = fig.add_subplot(111)
        else:
            ax = fig.axes[0]
        ax.plot(xax_, np.exp(log_dens), label=label)

        ax.fill_between(
            x=xax_,
            y1=np.exp(log_dens),
            where=(xax_ >= xax_[left_i + 1]) & (xax_ < xax_[right_i + 1]),
            # color="b",
            alpha=0.2,
        )
        ax.set_ylabel("PDF")
        ax.legend()

        # ax.axes.axvline(xax_[center_i], ls='--')

        # props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
        # place a text box in upper left in axes coords

        # ax.text(
        # 	0.05,
        # 	0.95,
        # 	textstr,
        # 	transform=ax.transAxes,
        # 	fontsize=14,
        # 	verticalalignment="top",
        # 	bbox=props,
        # )

    # test if area under the curve is "1"
    print(
        "Normalization test: ",
        quad(
            interp1d(xax_, np.exp(log_dens)),
            bounds[0],
            bounds[-1],
            epsabs=1e-5,
            epsrel=1e-5,
        ),
    )
    return xax_, log_dens, stats
