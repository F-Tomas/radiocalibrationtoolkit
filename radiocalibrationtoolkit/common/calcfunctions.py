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

NSIDE = 64
pixel_theta, pixel_phi = hp.pix2ang(NSIDE, np.arange(hp.nside2npix(NSIDE)))
ANGULAR_FACTOR = hp.pixelfunc.nside2pixarea(NSIDE) * np.sin(pixel_theta)
CONSTANT_FACTOR = (2 * constants.k / constants.c**2) * constants.physical_constants[
    "characteristic impedance of vacuum"
][0]


def voltage_squared_spectral_density(antenna_map: np.array, galactic_map: np.array, frequency_MHz: float) -> float:
    """
    Calculates the voltage squared spectral density.

    Args:
        antenna_map (np.array): Array of antenna data.
        galactic_map (np.array): Array of galactic data.
        frequency_MHz (float): Frequency in MHz.

    Returns:
        float: Voltage squared spectral density.
    """
    return (frequency_MHz * 1e6) ** 2 * CONSTANT_FACTOR * integrate_hpmap(galactic_map * antenna_map**2)


def add_integral_edges_to_x_axis(x: Union[list, np.ndarray], integral_edges: Union[list, np.ndarray]) -> np.ndarray:
    """
    Adds integral edges to the x-axis.

    Args:
        x: List or NumPy array of x-axis values.
        integral_edges: List or NumPy array of integral edges.

    Returns:
        NumPy array with the x-axis values including the integral edges.

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
    Integrates a function over a part of the array.

    Args:
        x: List or NumPy array of x-axis values.
        interp_func: Interpolating function.
        lower_integral_boundary: Lower boundary of the integration interval.
        upper_integral_boundary: Upper boundary of the integration interval.

    Returns:
        The value of the integral of the function over the specified interval.

    """
    value_filter = (x >= lower_integral_boundary) & (x <= upper_integral_boundary)
    return integrate.trapz(interp_func(x[value_filter]), x[value_filter])


def integrate_spectral_density(sd_DF, integrated_MHz_bands, integrating_method="quad", point=None):
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
        print("[ERROR] Not enough frequencies to integrate, calculate at least two frequency bins.")
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

    power_density_DF = pd.DataFrame(index=list(lst_range), columns=freq_Mhz_range, dtype="float64")

    temp = []
    for frequency_MHz in tqdm(freq_Mhz_range):

        galactic_map = hp.ma(hp.pixelfunc.ud_grade(galactic_map_inst.generate(frequency_MHz), NSIDE)).copy()
        # these to two aproaches are equivalent, but the second is ~25% faster
        # antenna_map_inst = antenna_inst.convert2hp(frequency=frequency_MHz, quantity='absolute', **update_antenna_conventions)
        antenna_map_inst = Data2hpmap(
            antenna_inst.get(frequency=frequency_MHz, quantity="absolute"), **update_antenna_conventions
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
                galactic_map = rotate_default_hpmap_array_from_galactic2local_coordinates(
                    galactic_map_in_galactic_cs,
                    lst,
                    latitude,
                    nside_transform=nside_transform,
                )
            else:
                galactic_map.mask = create_local_mask(NSIDE, rotation_parameters)
                antenna_map = antenna_map_inst.get_map(rotator=rotator)
            # todo fragmentation warn
            #  power_density_DF.loc[lst, frequency_MHz] =
            temp_row.append(
                voltage_squared_spectral_density(
                    antenna_map=antenna_map,
                    galactic_map=galactic_map,
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
    spectra_df = pd.DataFrame([i for i in time_trace_df.iloc[:, DSC:].apply(abs_fft, axis=1)])
    spectra_df.columns = np.arange(N // 2 + 1) / (N / 2) * sampling_frequency_MHz / 2
    return spectra_df


def get_bootstrap_CI(data, B=5000, alpha=1-0.6826, statistic_func=np.median):
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
