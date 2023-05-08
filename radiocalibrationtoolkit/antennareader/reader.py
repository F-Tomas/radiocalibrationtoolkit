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
    This class reads the antenna pattern from an XML file format and can interpolate the pattern for a finer spacing or convert it to a healpy format.

    Parameters
    ----------
    antennaPath : str
        The path of the XML file containing the antenna pattern data.
    frange : Tuple[float, float], optional
        The frequency range of interest, by default [0, np.inf].

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
    get(frequency: float, quantity: str, interp_phi: list = None, interp_theta: list = None, return_interp_func: bool = False) -> pd.DataFrame:
        Returns antenna pattern data as a pandas DataFrame for the given frequency and quantity. 
        Interpolates the data based on interp_phi and interp_theta if provided.
    
    Returns
    -------
    pd.DataFrame
        A pandas DataFrame of the antenna pattern data.
    """
    def __init__(self, antennaPath: str, frange:Tuple[float, float]=[0, np.inf]):
        self._xml_dict = collections.defaultdict(dict)
        root_ = GetXMLRoot(antennaPath)
        for topBranch in root_:
            # print(topBranch)
            if "idfreq" in topBranch.attrib.keys():
                # print('e')
                freq = float(topBranch.attrib["idfreq"])
                if (freq >= frange[0]) and (freq <= frange[1]):
                    self._xml_dict[topBranch.tag][freq] = convert_xml_string_to_floats(topBranch.text)
            else:
                values = np.asarray(list(dict.fromkeys(convert_xml_string_to_floats(topBranch.text))))
                if topBranch.tag == 'frequency':
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
        return_interp_func: bool = False
    ) -> pd.DataFrame:
        """
        Get the antenna pattern data for a given frequency and quantity.

        Parameters:
            frequency (float): The frequency for which to retrieve the antenna pattern data.
            quantity (str): The quantity to retrieve. Valid values are "absolute", "theta_amp", "phi_amp", "theta_phase", and "phi_phase".
            interp_phi (list, optional): A list of phi values for which to interpolate the data. Defaults to None.
            interp_theta (list, optional): A list of theta values for which to interpolate the data. Defaults to None.
            return_interp_func (bool, optional): Whether to return the interpolation function. Defaults to False.

        Returns:
            pd.DataFrame: A pandas dataframe containing the antenna pattern data for the chosen frequency and quantity. The rows are the phi angles and columns the theta angles.
        """
        if quantity not in self._quantity_list:
            self.show_quantities()
        if quantity == "absolute":
            # print(
            #    "[INFO] Trying to identify keys in the dictionary representing amplitudes of theta and phi polarization"
            # )
            # print("[INFO] This identification is based on which keys end with 'amp'")
            amp_keys = [k for k in self._xml_dict.keys() if k.endswith("amp")]
            if len(amp_keys) == 2:
                # print("[INFO] I found these keys: {} It should work.".format(amp_keys))
                values = np.sqrt(
                    (1 / 2)
                    * (
                        self._xml_dict[amp_keys[0]][frequency] ** 2
                        + self._xml_dict[amp_keys[1]][frequency] ** 2
                    )
                )
            else:
                print(
                    "[INFO] Tried to identify keys in the dictionary representing amplitudes of theta and phi polarization"
                )
                print(
                    "[INFO] This identification was based on which keys end with 'amp'"
                )
                print(
                    "[ERROR] No keys found, do your absolute value manually or fix the tags in the XML so that the polarization amplitudes end with 'amp'\
                      The formula is (1/2*theta_amp**2+phi_amp**2)^(1/2)  "
                )
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
        return temporary_DF

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
        if yes == True:
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
        add_invisible_sky: bool = True,
    ):
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
        # print(args)
        return AntennaPattern._Convert2hp(self, **args)

    class _Convert2hp(Data2hpmap):
        def __init__(self, AntennaPattern, **kwargs):
            self._frequency = kwargs["frequency"]
            self._df = AntennaPattern.get(
                frequency=self._frequency, quantity=kwargs["quantity"]
            )
            del_keys(kwargs, ["frequency", "quantity"])
            super().__init__(self._df, **kwargs)

        def get_frequency(self):
            return self._frequency

