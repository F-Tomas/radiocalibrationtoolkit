"""
This module contains three classes: GMOSS, ULSA, and SSM.

GMOSS - Generates the Galactic and extragalactic sky brightness temperature at frequencies ranging from 30 MHz to 80 MHz based on a table of precomputed values.

ULSA - Generates the ultra-large scale structure (ULSS) component of the sky brightness temperature at frequencies ranging from 30 MHz to 80 MHz based on a table of precomputed values.

SSM - Generates the sky temperature due to the Sunyaev-Zeldovich effect at a specified frequency and HEALPix resolution.

Dependencies: pandas, numpy, healpy, os, and functions from the common.py module.

Usage:
Instantiate a class and call its "generate" method with a frequency in MHz as an argument. The "generate" method will return an array of brightness temperature values.

Note: These classes assume that the user has the necessary data files and functions installed on their system.
"""
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# some classes, I will move this later to a separate py file

import pandas as pd
import numpy as np
import healpy as hp
import os
from pygdsm import (
    GlobalSkyModel2016,
    GlobalSkyModel,
    LowFrequencySkyModel,
    HaslamSkyModel,
)
from pylfmap import LFmap
from ..common import *

CMB_temp = 2.7255

try:
    import SSM.SSM as ssm_func
except ImportError as e:
    print("[ERROR] SSM is not installed!")
    pass


_path = os.path.dirname(os.path.abspath(__file__))


class GMOSS:
    """
    A class to generate the Galactic diffuse emission map from GMRT observations
    between 30-80 MHz based on the tabulated data.

    Attributes:
    -----------
    gmoss_map_DF : pandas.DataFrame
        Dataframe containing the GMRT observed spectra for Galactic diffuse emission
        between 30-80 MHz.

    """

    def __init__(self):
        self.gmoss_map_DF = pd.read_hdf(
            os.path.join(_path, "tabulated/gmoss_spectra_1-400MHz.h5")
        )
        self.gmoss_map_DF.columns = self.gmoss_map_DF.columns.astype(float) * 1000
        temp = [
            hp.reorder(self.gmoss_map_DF[col].values, n2r=True)
            for col in self.gmoss_map_DF.columns
        ]
        self.gmoss_map_DF.iloc[:, :] = np.asarray(temp).T

    def generate(self, frequency_MHz):
        """
        Generates the Galactic diffuse emission map for the given frequency.

        Parameters:
        -----------
        frequency_MHz : float
            The frequency value in MHz for which the Galactic diffuse emission map needs
            to be generated.

        Returns:
        --------
        numpy.ndarray
            The generated Galactic diffuse emission map for the given frequency.
        """
        if frequency_MHz in self.gmoss_map_DF.columns:
            return self.gmoss_map_DF[frequency_MHz].values + CMB_temp
        else:
            print(
                "[ERROR] Frequency {} was not tabulated! Run first GMOSS, update the table, and then call."
            )


class ULSA:
    """
    A class to generate the UTR-2 Low-Frequency Sky Atlas maps from the tabulated data
    based on the specified index type.

    Attributes:
    -----------
    ulsa_map_DF : pandas.DataFrame
        Dataframe containing the UTR-2 Low-Frequency Sky Atlas maps between 30-80 MHz.
    """

    def __init__(self, index_type: "str" = "freq_dependent_index"):
        index_types = [
            "freq_dependent_index",
            "constant_index",
            "direction_dependent_index",
        ]
        try:
            path = os.path.join(_path, "tabulated/ulsa_{}.h5".format(index_type))
            self.ulsa_map_DF = pd.read_hdf(path)
        except Exception:
            traceback.print_exc()
            print(
                "[ERROR] File {} not found. Did you use correct index type? Possible index types: {}".format(
                    path, index_types
                )
            )
        self.ulsa_map_DF.columns = self.ulsa_map_DF.columns.astype(float)

    def generate(self, frequency_MHz):
        """
        Parameters:
        -----------
        frequency_MHz : float
            The frequency value in MHz for which the UTR-2 Low-Frequency Sky Atlas map
            needs to be generated.

        Returns:
        --------
        numpy.ndarray
            The generated UTR-2 Low-Frequency Sky Atlas map for the given frequency.
        """
        if frequency_MHz in self.ulsa_map_DF.columns:
            return self.ulsa_map_DF[float(frequency_MHz)].values
        else:
            print(
                "[ERROR] Frequency {} was not tabulated! Run first ULSA, update the table, and then call."
            )


class SSM:
    """
    A class to generate the simulated sky maps from the Global Sky Model (GSM) based on
    the specified frequency and nside.
    """

    def __init__(self) -> None:
        pass

    def generate(self, frequency_MHz, nside=64):
        """
        Parameters:
        -----------
        frequency_MHz : float
            The frequency value in MHz for which the simulated sky map needs to be
            generated.

        nside : int, optional
            The value of nside to use for the generation of simulated sky map. The default
            value is 64.

        Returns:
        --------
        numpy.ndarray
            The generated simulated sky map from the Global Sky Model (GSM) for the given
            frequency and nside.
        """
        data = ssm_func(freq=frequency_MHz, nside=nside, verbose=False, name="del")
        os.remove("del.hdf5")
        return data[0]

def add_CMB_decorator(cls):
    """
    Decorator function to add the Cosmic Microwave Background (CMB) temperature to the output of a class.

    Parameters:
    -----------
    cls : class
        The original class to be decorated.

    Returns:
    --------
    class
        The decorated class that adds the CMB temperature to the output.

    Notes:
    ------
    The CMB temperature (2.7255 K) is added to the output generated by the decorated class.

    Example:
    --------
    @add_CMB_decorator
    class OriginalClass:
        def generate(self, freq):
            # Original implementation

    # Decorate the class
    ModifiedClass = add_CMB_decorator(OriginalClass)

    # Create an instance of the modified class
    modified = ModifiedClass()

    # Generate the output with added CMB temperature
    result = modified.generate(freq)

    """
    class AddCMB(cls):
        def generate(self, freq):
            return super().generate(freq) + CMB_temp
    return AddCMB

GlobalSkyModel = add_CMB_decorator(GlobalSkyModel)
HaslamSkyModel = add_CMB_decorator(HaslamSkyModel)
