"""
Module Description
------------------

The `antennareader` module provides functionalities for reading antenna pattern data
from an XML file format and performing various operations on the data. It includes a
class `AntennaPattern` that allows users to read antenna pattern data from an XML file,
interpolate the data for finer spacing, and convert it to a healpix format.

The module is designed to work with Pandas DataFrames and supports data manipulation,
interpolation, and conversion to healpix format, which is useful for antenna pattern
analysis and visualization.

It also offers methods to retrieve volume data for a given quantity and frequencies.
"""


#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import collections

import pandas as pd
import numpy as np
import plotly.express as px

from typing import Type, Tuple
from typing import Union
from collections import OrderedDict
from lxml import etree as ET
from io import StringIO
from scipy.interpolate import interp1d

from ..common import *


class AntennaPattern:
    """
    This class reads antenna pattern data from an XML file and provides methods to work
    with the pattern data, interpolate it, convert it to healpix format, and more.

    Parameters
    ----------
    antennaPath : str
        The path of the XML file containing the antenna pattern data.
    frange : Tuple[float, float], optional
        The frequency range of interest, by default [0, np.inf].
    new_antenna_conventions : dict, optional
        A dictionary of new antenna conventions to apply for coordinate adjustments,
        by default an empty dictionary.

    Attributes
    ----------
    _xml_dict : defaultdict
        A defaultdict containing the antenna pattern data parsed from the XML file.
    _quantity_list : list
        A list of the different quantities available in the antenna pattern data.
    show_quantities : method
        Prints the list of available quantities.

    Methods
    -------
    get_raw() -> dict:
        Returns the antenna pattern data as a dictionary.
    get(frequency: float, quantity: str, interp_phi: list = None, interp_theta: list = None,
        return_interp_func: bool = False, return_with_updated_conventions: bool = True) -> pd.DataFrame:
        Returns antenna pattern data as a pandas DataFrame for the given frequency and quantity.
        Interpolates the data based on interp_phi and interp_theta if provided.
    get_volumetric_dataset(quantity: str, frequencies: Union[List[float], np.ndarray, None] = None) -> pd.DataFrame:
        Get volume data for a given quantity and frequencies.
    convert2hp(frequency: float, shift_phi: float = 0, flip_theta: bool = False, flip_phi: bool = False, in_degrees: bool = True, quantity: str = "absolute", add_invisible_sky: bool = False) -> Data2hpmap:
        Convert the antenna pattern data to healpix format at a given frequency.

    Returns
    -------
    Data2hpmap
        A Data2hpmap object containing the antenna pattern data in healpix format.

    Notes
    -----
    This class provides functionality to work with antenna pattern data, including
    interpolation, coordinate adjustments, and conversion to healpix format.

    Example
    -------
    # Create an instance of AntennaPattern
    antenna_inst = AntennaPattern(antennaPath="path/to/antenna.xml")

    # Get antenna pattern data for a specific frequency and quantity
    data = antenna_inst.get(frequency=50, quantity="theta_amp")

    # Convert the antenna pattern data to healpix format
    healpix_data = antenna_inst.convert2hp(frequency=50)

    """

    def __init__(
        self,
        antennaPath: str,
        frange: Tuple[float, float] = (0, np.inf),
        new_antenna_conventions: dict = {},
    ):
        self._antenna_conventions = {
            "flip_theta": False,
            "flip_phi": False,
            "shift_phi": 0,
            "add_invisible_sky": False,
            "in_degrees": True,
        }
        self._default_settings_for_conventions = self._antenna_conventions.copy()

        updated_convention_keys = set(self._antenna_conventions) & set(
            new_antenna_conventions
        )
        self._antenna_conventions.update(new_antenna_conventions)
        if updated_convention_keys != set():
            self.__updated_conventions = True
            print(
                "[INFO] These conventions will be updated on output:",
                updated_convention_keys,
            )
        else:
            self.__updated_conventions = False

        """Constructor method for AntennaPattern class."""
        self._xml_dict = collections.defaultdict(dict)
        root_ = GetXMLRoot(antennaPath)
        for topBranch in root_:
            if "idfreq" in topBranch.attrib.keys():
                freq = float(topBranch.attrib["idfreq"])
                if (freq >= frange[0]) and (freq <= frange[1]):
                    self._xml_dict[topBranch.tag][freq] = convert_xml_string_to_floats(
                        topBranch.text
                    )
            else:
                values = np.asarray(
                    list(dict.fromkeys(convert_xml_string_to_floats(topBranch.text)))
                )
                if topBranch.tag == "frequency":
                    values = values[(values >= frange[0]) & (values <= frange[1])]
                self._xml_dict["parameters"][topBranch.tag] = values
        self._quantity_list = list(self._xml_dict.keys())
        self._quantity_list.remove("parameters")
        self._quantity_list = self._quantity_list + ["absolute"]
        self.show_quantities()

    def show_quantities(self):
        """
        Prints the available quantities of the antenna pattern data.
        """
        print(
            "[INFO] Your keys are: {} Use them as the 'quantity'".format(
                self._quantity_list
            )
        )

    def get_raw(self) -> dict:
        """
        Returns the antenna pattern data as a dictionary.

        Returns
        -------
        dict
            A dictionary containing the antenna pattern data.
        """
        return self._xml_dict

    def get(
        self,
        frequency: float,
        quantity: str,
        interp_phi: list = None,
        interp_theta: list = None,
        return_interp_func: bool = False,
        return_with_updated_conventions: bool = True,
    ) -> pd.DataFrame:
        """
        Get the antenna pattern data for a given frequency and quantity.

        Parameters
        ----------
        frequency : float
            The frequency for which to retrieve the antenna pattern data.
        quantity : str
            The quantity to retrieve. Valid values are "absolute", "theta_amp", "phi_amp", "theta_phase", and "phi_phase".
        interp_phi : list, optional
            A list of phi values for which to interpolate the data. Defaults to None.
        interp_theta : list, optional
            A list of theta values for which to interpolate the data. Defaults to None.
        return_interp_func : bool, optional
            Whether to return the interpolation function. Defaults to False.

        Returns
        -------
        pd.DataFrame
            A pandas DataFrame containing the antenna pattern data for the chosen frequency and quantity. The rows are the phi angles and columns are the theta angles.
        """
        if quantity not in self._quantity_list:
            self.show_quantities()
        if quantity == "absolute":
            amp_keys = [k for k in self._xml_dict.keys() if k.endswith("amp")]
            if len(amp_keys) == 2:
                values = np.sqrt(
                    (1 / 2)
                    * (
                        self._xml_dict[amp_keys[0]][frequency] ** 2
                        + self._xml_dict[amp_keys[1]][frequency] ** 2
                    )
                )
            else:
                print(
                    "[INFO] No keys found, do your absolute value manually or fix the tags in the XML so that the polarization amplitudes end with 'amp'"
                )
                print("[INFO] The formula is (1/2 * theta_amp**2 + phi_amp**2)^(1/2)")
        else:
            values = self._xml_dict[quantity][frequency]

        temporary_DF = pd.DataFrame(
            values.reshape(
                self._xml_dict["parameters"]["phi"].size,
                self._xml_dict["parameters"]["theta"].size,
            ),
            index=self._xml_dict["parameters"]["phi"],
            columns=self._xml_dict["parameters"]["theta"],
        )
        # if interpolate
        if (interp_phi is not None) or (interp_theta is not None):
            temporary_DF = self._interpolate(
                temporary_DF,
                interp_phi=interp_phi,
                interp_theta=interp_theta,
                return_interp_func=return_interp_func,
            )
        temporary_DF.index.name = "phi"
        temporary_DF.columns.name = "theta"

        if self.__updated_conventions and return_with_updated_conventions:
            return update_antenna_coordinates(temporary_DF, **self._antenna_conventions)
        else:
            return temporary_DF

    def get_volumetric_dataset(
        self, quantity: str, frequencies: Union[List[float], np.ndarray, None] = None
    ) -> pd.DataFrame:
        """
        Get volume data for a given quantity and frequencies.

        Parameters
        ----------
        quantity : str
            The quantity to retrieve volume data for.
        frequencies : Union[List[float], np.ndarray, None], optional
            Frequencies for which volume data is requested.
            It can be either a list or a numpy array of floats.
            If not provided (None), the method will use default frequencies from the `antenna_inst` object.

        Returns
        -------
        pd.DataFrame
            A Pandas DataFrame containing volume data for the specified quantity and frequencies.
            The DataFrame is indexed by the frequencies and may have multiple columns depending on the structure of the data retrieved.
        """
        if frequencies is None:
            frequencies = self._xml_dict["parameters"]["frequency"].astype(float)

        df_list = []
        for f in frequencies:
            df = self.get(frequency=f, quantity=quantity)
            df_list.append(df)

        volume_data_df = pd.concat(df_list, keys=frequencies)
        volume_data_df.index.names = ["freq", volume_data_df.index.names[1]]

        return volume_data_df

    def _interpolate(
        self,
        df,
        interp_phi=None,
        interp_theta=None,
        return_interp_func=False,
    ):
        pattern_interp2d_func = RegularGridInterpolator(
            (df.index.values, df.columns.values), df.values
        )
        if return_interp_func:
            return pattern_interp2d_func

        if interp_phi is None:
            interp_phi = df.index.values
        if interp_theta is None:
            interp_phi = df.columns.values

        PHI_new, THETA_new = np.meshgrid(interp_phi, interp_theta)

        new_values = pattern_interp2d_func((PHI_new.T.flatten(), THETA_new.T.flatten()))

        return pd.DataFrame(
            new_values.reshape(interp_phi.size, interp_theta.size),
            columns=interp_theta,
            index=interp_phi,
        )

    def _flip_array(self, arr, yes):
        if yes:
            return arr[::-1]
        else:
            return arr

    def convert2hp(
        self,
        frequency: float,
        shift_phi: float = 0,
        flip_theta: bool = False,
        flip_phi: bool = False,
        in_degrees: bool = True,
        quantity: str = "absolute",
        add_invisible_sky: bool = False,
    ) -> Data2hpmap:
        """
        Convert the antenna pattern data to healpix format at a given frequency.

        Parameters
        ----------
        frequency : float
            The frequency at which to extract the antenna pattern data.
        shift_phi : float, optional
            The amount of rotation around the z-axis, in degrees.
        flip_theta : bool, optional
            If True, the antenna pattern is flipped about the x-axis.
        flip_phi : bool, optional
            If True, the antenna pattern is flipped about the y-axis.
        in_degrees : bool, optional
            If True, the antenna pattern data is assumed to be in degrees.
            If False, the data is assumed to be in radians.
        quantity : str, optional
            The quantity to extract. Must be one of the keys of the antenna pattern data.
        add_invisible_sky : bool, optional
            If True, add the invisible sky pixels to the healpix map.
            If False, only the visible pixels are included in the map.

        Returns
        -------
        Data2hpmap
            A Data2hpmap object containing the antenna pattern data in healpix format.
        """
        args = locals().copy()
        del args["self"]
        return AntennaPattern._Convert2hp(self, **args)

    class _Convert2hp(Data2hpmap):
        def __init__(self, AntennaPattern, **kwargs):
            self._frequency = kwargs["frequency"]
            quantity = kwargs["quantity"]
            del_keys(kwargs, ["frequency", "quantity"])

            # if updated antenna settings are provided by the "convert2hp" method,
            # the ones from the class are not used,
            if kwargs == AntennaPattern._default_settings_for_conventions:
                # nothing in kwargs is provided (settings from the class constructor are used)
                print(
                    "[INFO] Using updated antenna conventions provided at the class constructor"
                )
                kwargs = self._antenna_conventions
            else:
                # some settings are provided, then use those, do not update the kwargs
                pass

            # here we must get the default df without updated conventions, since these are going to be updated
            # later in the Data2hpmap class
            self._df = AntennaPattern.get(
                frequency=self._frequency,
                quantity=quantity,
                return_with_updated_conventions=False,
            )
            super().__init__(self._df, **kwargs)

        def get_frequency(self):
            return self._frequency
