"""
Module Description:
--------------------

The `mocktracegenerator` module provides functionalities for generating mock voltage
time traces. It includes a class `Mock_trace_generator` that allows users to generate
mock voltage traces based on voltage-to-density conversion values as a function of
sidereal time and frequency.

The module uses interpolation functions for the hardware response and impedance, along
with additional noise, to create realistic mock time traces. It also supports adding
non-background interference (NBI) to the traces and temperature-varying additional noise.

Users can specify various parameters, such as the size of the time trace, sampling
frequency, temperature range, and more, to customize the generated mock traces. The
module is useful for testing and simulating voltage time traces for antenna pattern
analysis and other related applications.
"""
#!/usr/bin/env python
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
from typing import Union, List, Tuple
from scipy.interpolate import RegularGridInterpolator
from scipy.constants import constants
from scipy.fft import ifft, fft
from scipy.signal import savgol_filter
import scipy.stats as stats
from ..common import *


class Mock_trace_generator:
    """
    A class that generates mock voltage time traces.

    Parameters
    ----------
    sidereal_voltage2_density_DF : pandas.DataFrame
        Dataframe containing the voltage-to-density conversion values as a function of sidereal time and frequency.
    hw_response_func : interp1d
        Interpolator function for the hardware response as a function of frequency.
    impedance_func : interp1d
        Interpolator function for the impedance as a function of frequency.
    voltage2ADC : int, optional
        The maximum voltage (in volts) corresponds to half of the ADC range. Default is 2048.
    time_trace_size : int, optional
        The size of the time trace. Default is 2048.
    sampling_frequency_MHz : float, optional
        The sampling frequency of the trace in MHz. Default is 250.

    Returns
    -------
    None
    """
    def __init__(
        self,
        sidereal_voltage2_density_DF: pd.DataFrame,
        hw_response_func: interp1d,
        impedance_func: interp1d,
        voltage2ADC: int = 2048,
        time_trace_size: int = 2048,
        sampling_frequency_MHz: float = 250,
    ) -> None:

        # read in voltage square density in antenna, not in HW!
        self._voltage2_density_interp2d_func = RegularGridInterpolator(
            (
                sidereal_voltage2_density_DF.index.values,
                sidereal_voltage2_density_DF.columns.values,
            ),
            sidereal_voltage2_density_DF.values,
            fill_value=None,
            bounds_error=False,
        )

        self._freq_MHz_bins = (
            np.arange(time_trace_size // 2 + 1)
            / (time_trace_size // 2)
            * (sampling_frequency_MHz / 2)
        )
        self._hw_response = hw_response_func(self._freq_MHz_bins)
        self._impedance = impedance_func(self._freq_MHz_bins)
        self._time_trace_size = time_trace_size
        self._sampling_frequency_MHz = sampling_frequency_MHz
        self._voltage2ADC = voltage2ADC

    def get_frequency_bins(self):
        """
        Return the frequency bins used in generating the mock trace.

        Returns
        -------
        numpy.ndarray
            An array of frequency values in MHz.
        """
        return self._freq_MHz_bins

    def generate_mock_trace(self, size:int=1, return_debug_dict:bool=False, *args, **kwargs) -> Union[pd.DataFrame, dict]:
        """
        Generate a mock voltage time trace.

        Parameters
        ----------
        size : int, optional
            The number of traces to generate. Default is 1.
        return_debug_dict : bool, optional
            If True, returns a dictionary with debug information in addition to the time trace.
            Default is False.
        '*'args : any
            Additional positional arguments to be passed to the '_generate_mock_trace' method.
        '**'kwargs : any
            Additional keyword arguments to be passed to the '_generate_mock_trace' method.

        Returns
        -------
        pd.DataFrame or dict
            If return_debug_dict is False, returns a pandas DataFrame containing the time trace data.
            If return_debug_dict is True, returns a dictionary containing both the time trace and debug information.
        """
        outputs = [] if return_debug_dict else pd.DataFrame()
        for i in tqdm(range(size)):
            out = self._generate_mock_trace(return_debug_dict=return_debug_dict, *args, **kwargs)
            if return_debug_dict:
                outputs.append(out)
            else:
                outputs = pd.concat((outputs, (out)))
        if isinstance(outputs, pd.DataFrame):
            outputs = outputs.reset_index(drop=True)
        return outputs

    def _generate_mock_trace(
        self,
        lst: Union[float, List[Tuple[int, int]]] = [0, 24],
        temp_celsius: Union[float, Tuple[float, float]] = [-10, 50],
        additional_noise: Union[float, interp1d] = 0,
        nbi: dict = {},
        nbi_err: float = 0,
        return_debug_dict: bool = False,
        rounding: bool = True,
        temperature_varying_additional_noise: bool = False
    ) -> pd.DataFrame:

        if isinstance(lst, list):
            lst_ = random.uniform(lst[0], lst[1])
        else:
            lst_ = lst

        if isinstance(temp_celsius, list):
            temp_K = random.uniform(temp_celsius[0], temp_celsius[1]) + 273.15
        else:
            temp_K = temp_celsius + 273.15

        freq_MHz_bins = self._freq_MHz_bins
        sidereal_voltage2_density = self._voltage2_density_interp2d_func(
            ([lst_], freq_MHz_bins)
        )

        thermal_noise_voltage2_density = 4 * constants.Boltzmann * temp_K * self._impedance

        total_voltage2_density = sidereal_voltage2_density + thermal_noise_voltage2_density

        total_voltage2_density_scattered = np.asarray([np.random.exponential(mu, size=1)[0] for mu in total_voltage2_density])
        
        total_voltage2_density_scattered_in_HW = (
            total_voltage2_density_scattered * self._hw_response
        )
        
        if isinstance(additional_noise, interp1d):
            additional_noise = additional_noise(self._freq_MHz_bins)
            # to prevent negative values from imperfect extrapolation
            additional_noise[additional_noise < 0] = 0
        if temperature_varying_additional_noise:
            additional_noise *= temp_K

        total_voltage2_density_scattered_in_HW_with_additional_noise = total_voltage2_density_scattered_in_HW + np.random.exponential(additional_noise, size=freq_MHz_bins.size)
        
        spectrum_voltage_density_in_HW_with_additional_noise = np.sqrt(total_voltage2_density_scattered_in_HW_with_additional_noise*self._time_trace_size * self._sampling_frequency_MHz * 1e6 / 2)

        # add NBI
        spectrum_voltage_density_in_HW_with_additional_noise_with_NBI = spectrum_voltage_density_in_HW_with_additional_noise.copy()
        
        for key in nbi:
            if float(key) > np.min(freq_MHz_bins) and float(key) < np.max(
                freq_MHz_bins
            ):
                spectrum_voltage_density_in_HW_with_additional_noise_with_NBI[
                    find_closest_index(freq_MHz_bins, float(key))
                ] = add_deviation(nbi[key], nbi_err)
            else:
                print("[WARN] NBI frequency outside the band!")

        spectrum_voltage_density_in_HW_with_additional_noise_with_NBI_in_complex = one_sided_2_complex_two_sided(
            spectrum_voltage_density_in_HW_with_additional_noise_with_NBI
        )

        mock_time_trace = np.real(ifft(spectrum_voltage_density_in_HW_with_additional_noise_with_NBI_in_complex)) * self._voltage2ADC + 0
        if rounding:
            mock_time_trace = np.asarray([round(i) for i in mock_time_trace])

        if return_debug_dict:
            return {
                "freq_MHz_bins": freq_MHz_bins,
                "sidereal_voltage2_density": sidereal_voltage2_density,
                "thermal_noise_voltage2_density": thermal_noise_voltage2_density,
                "total_voltage2_density_scattered": total_voltage2_density_scattered,
                'total_voltage2_density_scattered_in_HW':total_voltage2_density_scattered_in_HW,
                'total_voltage2_density_scattered_in_HW_with_additional_noise':total_voltage2_density_scattered_in_HW_with_additional_noise,
                "spectrum_voltage_density_in_HW": np.sqrt(total_voltage2_density_scattered_in_HW*self._time_trace_size * self._sampling_frequency_MHz * 1e6 / 2),
                "spectrum_voltage_density_in_HW_with_additional_noise": spectrum_voltage_density_in_HW_with_additional_noise,
                "spectrum_voltage_density_in_HW_with_additional_noise_with_NBI": spectrum_voltage_density_in_HW_with_additional_noise_with_NBI,
                "time_trace":mock_time_trace,
                "additional_noise":np.random.exponential(additional_noise, size=freq_MHz_bins.size)
            }
        else:
            return pd.DataFrame(
                [[lst_, temp_K - 273.15, *mock_time_trace]],
                columns=["lst", "temp_c", *np.asarray(range(self._time_trace_size))/self._sampling_frequency_MHz],
            )

