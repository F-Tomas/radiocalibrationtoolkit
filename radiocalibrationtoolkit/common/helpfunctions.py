#!/usr/bin/env python
# -*- coding: utf-8 -*-


import random
import numpy as np
import healpy as hp
import pandas as pd
import re
import os
import collections
import argparse
from pathlib import Path
from typing import Type, Tuple, Optional, Union, List
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import quad
from scipy.fft import fft, ifft
from sklearn.neighbors import KernelDensity
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.axes as axes
import matplotlib.pyplot as plt
from healpy.rotator import Rotator
from healpy.newvisufunc import projview

from collections import OrderedDict
from lxml import etree as ET
from io import StringIO


def compare_maps(
    map_dict: dict,
    main_title: str = "",
    rotation_parameters: list = None,
    coord: list = ["G"],
    mask: list = False,
    projection_type: str = "mollweide",
    cmap: str = "jet",
    show_plot_comparison: bool = True,
    verbose: bool = False,
    figsize: Tuple = (6, 6.5),
) -> dict:
    """
    Compare maps in a dictionary using matplotlib and return an integrated comparison dictionary.

    Parameters:
    -----------
    map_dict : dict
        A dictionary containing maps to compare.

    main_title : str, optional (default="")
        The title for the entire comparison plot.

    rotation_parameters : list, optional (default=None)
        List of longitude and latitude values for rotating the plot.

    coord : list, optional (default=['G'])
        Coordinate system to be used. Can be one of ['G', 'H', 'C', 'S'].

    mask : list, optional (default=False)
        Whether to apply a mask to the maps.

    projection_type : str, optional (default="mollweide")
        The type of projection to use for the plot.

    cmap : str, optional (default="jet")
        The color map to use for the plot.

    show_plot_comparison : bool, optional (default=True)
        Whether to display the plot.

    verbose : bool, optional (default=False)
        Whether to print verbose output during the function execution.

    Returns:
    --------
    dict
        An integrated comparison dictionary.

    Raises:
    -------
    None

    """
    size = len(map_dict.keys())
    keys = list(map_dict.keys())
    i = 1
    integrated_comparison_dict = {}
    if show_plot_comparison:
        fig = plt.figure(figsize=figsize)
    for row in range(size):
        for col in range(size):
            if verbose:
                print(i, row, col)
            if row == col:
                m = map_dict[keys[row]]
                m.mask = mask
                norm = "log"
                title = keys[row]
                min_t = 3500
                max_t = 35000
                xlabel = ""
                cbar = False
            else:
                m_row = map_dict[keys[row]]
                m_col = map_dict[keys[col]]
                m_row.mask = mask
                m_col.mask = mask
                m = m_row / m_col - 1
                norm = None
                title = keys[row] + "/" + keys[col]
                min_t = -0.5
                max_t = 0.5
                integrated_ratio = integrate_hpmap(m_row) / integrate_hpmap(m_col) - 1
                xlabel = f"{integrated_ratio:+.2f}"
                cbar = True
            if col > row:
                integrated_comparison_dict[title] = [integrated_ratio]
            if col == 0:
                ylabel = keys[row]
            else:
                ylabel = ""
            if row == 0:
                title = keys[col]
            else:
                title = ""
            if show_plot_comparison:
                projview(
                    m,
                    coord=coord,
                    rot=rotation_parameters,
                    norm=norm,
                    graticule=False,
                    graticule_labels=False,
                    unit="",
                    cb_orientation="vertical",
                    min=min_t,
                    max=max_t,
                    # latitude_grid_spacing=30,
                    projection_type=projection_type,
                    title=title,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    xtick_label_color="white",
                    cmap=cmap,
                    sub=(size, size, i),
                    fontsize={"xlabel": 10, "ylabel": 12, "title": 12},
                    fig=fig,
                    cbar=False
                    # cbar=cbar,
                )
                ax = plt.gca()
                ax.set_title(title, size=12, weight="bold")
                if i == 2:
                    cax = fig.add_axes(
                        [0.2, 0.0, 0.5, 0.025]
                    )  # [left, bottom, width, height]
                    cbar = plt.colorbar(
                        # label="Colorbar Label",
                        cmap=cmap,
                        norm=plt.Normalize(vmin=-0.5, vmax=0.5),
                        cax=cax,
                        orientation="horizontal",
                        extend="both",
                    )
                    cbar.set_ticks(
                        [min_t, 0, max_t]
                    )  # Replace with desired tick values
            else:
                pass
            i += 1
    if show_plot_comparison:
        plt.suptitle(main_title, fontsize=14)
        plt.subplots_adjust(top=0.9, wspace=0.05, hspace=0.1)
    return integrated_comparison_dict


def create_local_mask(
    nside: int, rotation_parameters: Tuple[float, float]
) -> np.ndarray:
    """
    Creates a mask that is True for pixels above the equator and False for those below.

    Parameters
    ----------
    nside : int
        The nside parameter for the Healpix map.
    rotation_parameters : tuple of float
        A tuple containing the Galactic longitude and latitude (in degrees) used to rotate the map.

    Returns
    -------
    np.ndarray
        A boolean array with True for pixels above the equator and False for those below.
    """
    m = hp.ma(np.arange(hp.nside2npix(nside), dtype=np.double))
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    pixel_theta, pixel_phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))

    r = Rotator(coord=["G", "C"], rot=rotation_parameters)
    pixel_theta, pixel_phi = r(pixel_theta, pixel_phi)

    mask[pixel_theta > np.pi / 2] = 1
    return mask


def print_map_properties(map_: np.ndarray) -> None:
    """
    Prints properties of a Healpix map.

    Parameters
    ----------
    map : np.ndarray
        A one-dimensional array representing the map.

    Returns
    -------
    None
    """
    print(
        map_.size,
        hp.get_nside(map_),
        np.rad2deg(hp.pixelfunc.nside2resol(hp.get_nside(map_))),
    )


def find_index_on_CDF(log_dens: np.ndarray, X: np.ndarray, p: float) -> int:
    """
    Returns the index of the array `X` such that the cumulative distribution function (CDF) of `log_dens`
    evaluated on `X` is closest to `p`.

    Parameters
    ----------
    log_dens : np.ndarray
        Array of log-density values.
    X : np.ndarray
        The array of values to evaluate the CDF on.
    p : float
        The probability threshold.

    Returns
    -------
    int
        The index of the array `X` such that the CDF of `log_dens` evaluated on `X` is closest to `p`.
    """
    return np.argmin(np.abs(np.cumsum(np.exp(log_dens) * np.diff(X.flatten())[0]) - p))


def create_rotation_parameters(lst: float, latitude: float) -> Tuple[float, float]:
    """
    Creates a tuple of Galactic longitude and latitude (in degrees) from local sidereal time and latitude.

    Parameters
    ----------
    lst : float
        The local sidereal time in hours.
    latitude : float
        The latitude of the observer in degrees.

    Returns
    -------
    tuple of float
        A tuple containing the Galactic longitude and latitude (in degrees).
    """
    return (180 + (lst * 15)) % 360, -(latitude - 90)


# def create_rotator(lst, latitude, coord=["G"], inv=False):
def create_rotator(
    lst: Optional[float],
    latitude: Optional[float],
    coord: List[str] = ["G"],
    inv: bool = False,
) -> Rotator:
    """
    Creates a Rotator object for converting between Galactic and celestial (e.g., equatorial) coordinates.

    Parameters
    ----------
    lst : float or None
        The local sidereal time in hours. If None, defaults to 0.
    latitude : float or None
        The latitude of the observer in degrees. If None, defaults to 0.
    coord : list of str, optional
        A list of two strings indicating the input and output coordinate systems. Defaults to ["G"] (Galactic).
    inv : bool, optional
        Whether to perform an inverse rotation. Defaults to False.

    Returns
    -------
    Rotator
        A Rotator object that can be used to convert between Galactic and celestial coordinates.
    """
    if (lst == None) or (latitude == None):
        return Rotator(coord=coord)
    else:
        return Rotator(
            coord=coord, rot=create_rotation_parameters(lst, latitude), inv=inv
        )


# def rotate_default_hpmap_array_from_galactic2local_coordinates(m, lst, latitude, inv=False, nside_transform = 512):
def rotate_default_hpmap_array_from_galactic2local_coordinates(
    m: np.ndarray,
    lst: Optional[float],
    latitude: Optional[float],
    inv: bool = False,
    nside_transform: int = 512,
) -> np.ndarray:
    """
    Rotates a Healpix map from Galactic to local coordinates.

    Parameters
    ----------
    m : np.ndarray
        A one-dimensional array representing the map.
    lst : float or None
        The local sidereal time in hours. If None, defaults to 0.
    latitude : float or None
        The latitude of the observer in degrees. If None, defaults to 0.

    """
    nside_original = hp.npix2nside(len(m))
    m = hp.ma(hp.pixelfunc.ud_grade(m, nside_transform))
    pixel_theta, pixel_phi = hp.pix2ang(
        nside_transform, np.arange(hp.nside2npix(nside_transform))
    )
    rotator = create_rotator(lst, latitude, coord=["G", "C"], inv=bool(1 - inv))
    pixel_rot_theta, pixel_rot_phi = rotator(pixel_theta, pixel_phi)
    new_pix_order = hp.ang2pix(nside_transform, pixel_rot_theta, pixel_rot_phi)
    m = m[new_pix_order]
    m = hp.ma(hp.pixelfunc.ud_grade(m, nside_original))
    return m


def my_modulo(arr, modulo):
    """
    Compute the modulo of each element in an array, except for elements that are equal to the modulo.

    Parameters
    ----------
    arr : array_like
        The input array.
    modulo : scalar
        The modulo value.

    Returns
    -------
    ndarray
        The output array with the same shape as the input array.

    Examples
    --------
    >>> my_modulo([1, 2, 3, 4, 5, 6, 7], 3)
    array([1, 2, 0, 1, 2, 0, 1])
    """
    return np.asarray([x if x == modulo else x % modulo for x in arr])


def mkdir(directory):
    """
    Create a directory if it doesn't already exist.

    Parameters:
    -----------
    directory : str
        Directory path to be created.

    Returns:
    --------
    None

    Notes:
    ------
    - This function checks if the directory already exists.
    - If the directory doesn't exist, it creates a new directory.
    - If the directory already exists, it prints a message indicating that the directory already exists.

    Example:
    --------
    # Example usage
    mkdir("path/to/directory")
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
        print("Directory created successfully")
    else:
        print("Directory already exists")


class Data2hpmap:
    """
    Interpolates a pandas DataFrame onto a healpix map.
    Supports half-sky and full-sky DataFrames, and includes options to shift and flip the data.

    Parameters
    ----------
    data : pandas DataFrame or 2D array
        The input data. If a pandas DataFrame is provided, the index and column names will be used
        as the phi and theta values, respectively. If a 2D array is provided, `phi` and `theta`
        must be specified.
    phi : array_like, optional
        The phi values corresponding to the columns of the input data. Ignored if `data` is a
        pandas DataFrame. Default is None.
    theta : array_like, optional
        The theta values corresponding to the rows of the input data. Ignored if `data` is a
        pandas DataFrame. Default is None.
    shift_phi : float, optional
        The amount to shift the phi values by (in radians). Default is 0.
    flip_theta : bool, optional
        Whether to flip the theta values. Default is False.
    flip_phi : bool, optional
        Whether to flip the phi values. Default is False.
    in_degrees : bool, optional
        Whether the input phi and theta values are in degrees. Default is False.
    add_invisible_sky : bool, optional
        Whether to add an invisible sky below the input data (for half-sky maps). Default is False.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If the input data is not a pandas DataFrame and `phi` and `theta` are not provided.

    Notes
    -----
    This class uses the `scipy.interpolate.RegularGridInterpolator` function to interpolate the data onto
    a healpix map. The interpolated map can be obtained using the `get_map` method.

    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> data = np.random.randn(10, 20)
    >>> phi = np.linspace(0, 2*np.pi, 20, endpoint=False)
    >>> theta = np.linspace(0, np.pi/2, 10)
    >>> df = pd.DataFrame(data, index=theta, columns=phi)
    >>> hp_map = Data2hpmap(df, shift_phi=np.pi/4).get_map(nside=64)
    """

    def __init__(
        self,
        data,
        phi=None,
        theta=None,
        shift_phi=0,
        flip_theta=False,
        flip_phi=False,
        in_degrees=False,
        add_invisible_sky=False,
    ):
        if isinstance(data, pd.DataFrame):
            df = data.copy(deep=True)
        else:
            if isinstance(phi, (list, np.ndarray)) and isinstance(
                theta, (list, np.ndarray)
            ):
                df = pd.DataFrame(data, columns=theta, index=phi)
            else:
                print("[ERROR] Data are not a pandas dataframe. Define theta= and phi=")

        # just in case...
        df.columns = df.columns.values.astype(float)
        df.index = df.index.values.astype(float)

        if in_degrees:
            df.columns = np.deg2rad(df.columns.values)
            df.index = np.deg2rad(df.index.values)
            shift_phi = np.deg2rad(shift_phi)

        # reorder to ascending order if needed
        if df.columns[0] > df.columns[-1]:
            df = df.loc[:, ::-1]
        if df.index[0] > df.index[-1]:
            df = df.loc[::-1, :]

        if add_invisible_sky:
            if flip_theta:
                df.iloc[:, :] = df.iloc[:, ::-1].values
            if df.columns[0] > np.pi / 2:  # in case that theta goes from pi/2 to pi
                df.columns = df.columns - np.pi / 2
            elif df.columns[0] < 0:
                print("[ERROR] Theta for the half sky cannot be negative.")
                return None
            invisible_DF = pd.DataFrame(
                np.zeros(df.shape)[:, :-1],
                columns=-df.columns.values[1:][::-1],
                index=df.index,
            )
            df_complete = pd.concat((invisible_DF, df), axis=1)
        else:
            df_complete = df
            if flip_theta:
                df_complete.iloc[:, :] = df_complete.iloc[:, ::-1].values

        phi_g = df_complete.index.values
        theta_g = df_complete.columns.values

        # adjust to conventions
        if phi_g.max() > np.pi:  # if phi range is 0,2pi change it to -/+pi
            phi_g = phi_g - np.pi
        if theta_g.min() < 0:  # if theta is -/+pi/2 change it to 0,pi
            theta_g = theta_g + np.pi / 2

        # shifting in azimuth
        if shift_phi != 0:
            shifted = my_modulo((phi_g + np.pi + shift_phi), (2 * np.pi))
            new_row_order = np.argsort(np.argsort(shifted))
            df_complete.iloc[:, :] = df_complete.iloc[new_row_order, :].values

        if flip_phi:
            df_complete.iloc[:, :] = df_complete.iloc[::-1, :].values

        # return df_complete
        # internally healpy has different ordering, the 'theta' is in descending order, starting with pi, flipping it here
        self._interp2d_func = RegularGridInterpolator(
            (phi_g, theta_g), df_complete.values[:, ::-1]
        )
        self.df_complete = df_complete

    def get_grid_df(self) -> pd.DataFrame:
        """
        Returns the input data as a pandas DataFrame, with the rows corresponding to the `phi` values and
        the columns corresponding to the `theta` values.

        Returns
        -------
        pd.DataFrame
            The input data as a pandas DataFrame.
        """
        return self.df_complete

    def get_map(
        self, nside: int = 64, rotator: type(Rotator) = Rotator(rot=[0, 0])
    ) -> np.ndarray:
        """
        Interpolate the input data onto a full-sky map with the specified `nside` and `rotator`.

        Parameters
        ----------
        nside : int, optional
            The resolution parameter for the HEALPix map. Defaults to 64.
        rotator : type(Rotator), optional
            The rotator to use to rotate the interpolated map. Defaults to Rotator(rot=[0, 0]).

        Returns
        -------
        np.ndarray
            The interpolated HEALPix map.
        """
        # find theta and phi that need to be interpolated
        pixel_theta, pixel_phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
        pixel_rot_theta, pixel_rot_phi = rotator(pixel_theta, pixel_phi)
        # print(np.min(pixel_rot_theta), np.max(pixel_rot_theta), np.min(pixel_rot_phi), np.max(pixel_rot_phi) )
        try:
            return self._interp2d_func((pixel_rot_phi, pixel_rot_theta))
        except Exception:
            traceback.print_exc()
            print(
                "[ERROR] Is the input map for the whole sky? If it is only for half, use 'add_invisible_sky=True'"
            )

    def get_interp_func(self):
        """
        Get the interpolation function.

        Returns
        -------
        callable
            The interpolation function.

        Notes
        -----
        The interpolation function is a `RegularGridInterpolator` object that can be used to interpolate the data onto a
        regular grid. The function takes a tuple (phi, theta) as input, where `phi` and `theta` are arrays of the same shape,
        and returns the interpolated values at the corresponding (phi, theta) positions.
        """
        return self._interp2d_func


def fig_to_image(fig):
    """
    Converts a Matplotlib figure object to a NumPy array of RGB values.

    Parameters:
    -----------
    fig : matplotlib.figure.Figure object
        Matplotlib figure object to be converted.

    Returns:
    --------
    image_from_plot : ndarray
        3D NumPy array of RGB values representing the image.
    """
    canvas = FigureCanvas(fig)
    ax = fig.gca()
    canvas.draw()  # draw the canvas, cache the renderer

    image_from_plot = np.frombuffer(canvas.tostring_rgb(), dtype="uint8")
    return image_from_plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))


def convert_xml_string_to_floats(xml_string):
    """
    Converts an XML string containing whitespace-separated float values to a NumPy array of float values.

    Parameters:
    -----------
    xml_string : str
        XML string containing whitespace-separated float values.

    Returns:
    --------
    floats : ndarray
        NumPy array of float values extracted from the XML string.
    """
    temp = re.split("\s+", xml_string)
    return np.asarray([v for v in temp if v != ""]).astype(float)


def GetXMLRoot(XMLFile):
    """
    Reads an XML file and returns the root element of the parsed XML tree.

    Parameters:
    -----------
    XMLFile : str
        Path to the XML file to be read.

    Returns:
    --------
    root : xml.etree.ElementTree.Element
        Root element of the parsed XML tree.
    """
    fstring = [fstring for fstring in open(XMLFile, "r").readlines() if fstring.strip()]
    fstring = "".join(fstring)

    # to create proper root
    fstring = "<rp>" + fstring + "</rp>"
    parser = ET.XMLParser(remove_blank_text=True, remove_comments=True, recover=True)

    root = ET.fromstring(fstring, parser)
    root.tag
    root.attrib

    return root


def read_hw_file(
    hw_file_path: str,
    return_tabulated: bool = False,
    return_interp_func: bool = True,
    interp_args: dict = {},
) -> dict:
    """
    Reads an XML file containing hierarchical data, extracts data into a nested dictionary of NumPy arrays and/or interpolation functions.

    Parameters:
    -----------
    hw_file_path : str
        Path to the XML file to be read.
    return_tabulated : bool, optional
        If True, returns data in a nested dictionary of NumPy arrays, else returns data in a nested dictionary of interpolation functions.
        Default is False.
    return_interp_func : bool, optional
        If True, returns data in a nested dictionary of interpolation functions, else returns data in a nested dictionary of NumPy arrays.
        Default is True.
    interp_args : dict, optional
        Dictionary of optional arguments to be passed to the `interp1d` function of NumPy.
        Default is an empty dictionary.

    Returns:
    --------
    xml_dict : dict
        A nested dictionary containing the extracted data. If `return_tabulated` is True, values in the dictionary are NumPy arrays.
        Else, values are interpolation functions.
    """
    if return_tabulated == True:
        return_interp_func = False
    nested_dict = lambda: collections.defaultdict(nested_dict)
    xml_dict = nested_dict()
    root_ = GetXMLRoot(hw_file_path)
    for top_branch in root_:
        print(top_branch)
        for sub_branch in top_branch:
            for ssub_branch in sub_branch:
                if return_interp_func:
                    x = convert_xml_string_to_floats(ssub_branch.getchildren()[0].text)
                    y = convert_xml_string_to_floats(ssub_branch.getchildren()[1].text)
                    xml_dict[sub_branch.tag][sub_branch.attrib["id"]] = interp1d(
                        x, y, **interp_args
                    )
                if return_tabulated:
                    for sssub_branch in ssub_branch:
                        xml_dict[sub_branch.tag][sub_branch.attrib["id"]][
                            sssub_branch.tag
                        ] = convert_xml_string_to_floats(sssub_branch.text)
    return xml_dict


def hpmap2grid(m: List, xsize: int = 1000) -> Tuple:
    """
    Interpolate a HEALPix map onto a regular grid.

    Parameters:
    -----------
    m : List
        HEALPix map array.
    xsize : int, optional (default=1000)
        The size of the x-axis of the interpolated map.

    Returns:
    --------
    PHI : numpy.ndarray
        2D array of azimuthal angles (longitude).
    THETA : numpy.ndarray
        2D array of polar angles (colatitude).
    grid_map : numpy.ndarray
        Interpolated map.
    """
    ysize = xsize // 2
    theta = np.linspace(np.pi, 0, ysize)
    phi = np.linspace(-np.pi, np.pi, xsize)

    PHI, THETA = np.meshgrid(phi, theta)

    # THETA, PHI = r(THETA.flatten(), PHI.flatten())
    THETA = THETA.reshape(ysize, xsize)
    PHI = PHI.reshape(ysize, xsize)

    nside = hp.pixelfunc.npix2nside(len(m))
    grid_pix = hp.ang2pix(nside, THETA, PHI, nest=False)
    grid_map = m[grid_pix]
    return PHI, THETA, grid_map[::-1, ::-1]


def integrate_hpmap(m: List) -> float:
    """
    Integrate the HEALPix map over the sky.

    Parameters:
    -----------
    m : List
        HEALPix map array.

    Returns:
    --------
    flux : float
        The integrated flux of the map.
    """
    PHI, THETA, grid_map = hpmap2grid(m)
    phi_delta = abs(np.diff(PHI[0, :])[0])
    phi_theta = abs(np.diff(THETA[:, 0])[0])
    return np.sum(phi_delta * phi_theta * grid_map * np.sin(THETA))


def del_keys(d: dict, keys: Union[str, List[str]]) -> None:
    """
    Remove one or more keys from a dictionary.

    Parameters:
    -----------
    d : dict
        The dictionary to modify.
    keys : str or List[str]
        The key(s) to remove from the dictionary.
    """
    if isinstance(keys, list) != True:
        keys = [keys]
    for key in keys:
        del d[key]


def add_deviation(mu: float, error: float) -> float:
    """
    Add a random deviation to a value based on a fractional error.

    Parameters:
    -----------
    mu : float
        The mean value.
    error : float
        The fractional error of the value.

    Returns:
    --------
    new_mu : float
        The new mean value with the deviation added.
    """
    std = error * mu
    return np.random.normal(mu, std)


def one_sided_2_two_sided(spectrum: np.ndarray, sign: int = 1) -> np.ndarray:
    """
    Convert a one-sided power spectrum into a two-sided power spectrum.

    Parameters:
    -----------
    spectrum : numpy.ndarray
        One-sided power spectrum.
    sign : int, optional (default=1)
        The sign to use when computing the negative frequencies in the two-sided spectrum.

    Returns:
    --------
    two_sided_spectrum : numpy.ndarray
        Two-sided power spectrum.
    """
    return np.hstack((spectrum, sign * spectrum[1:][::-1][1:]))


def generate_one_sided_random_phase(size: int) -> np.ndarray:
    """
    Generate one-sided random phase for a given size.

    Parameters
    ----------
    size : int
        Size of the phase array.

    Returns
    -------
    np.ndarray
        One-sided random phase array.

    """
    random_phase = np.random.uniform(low=np.pi, high=-np.pi, size=size)
    #  phase_ramp = np.linspace(np.pi, -np.pi, size)
    phase_ramp = np.linspace(-np.pi, np.pi, size)
    random_ramped_phase = random_phase + phase_ramp * 0
    random_ramped_phase_wrapped = np.unwrap(random_ramped_phase)
    random_ramped_phase_wrapped -= (
        2 * np.pi * np.floor((random_ramped_phase_wrapped + np.pi) / (2 * np.pi))
    )
    return random_ramped_phase_wrapped


def one_sided_2_complex_two_sided(
    one_sided_abs: np.ndarray, one_sided_phase: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Convert one-sided amplitude and phase to complex two-sided array.

    Parameters
    ----------
    one_sided_abs : np.ndarray
        One-sided amplitude array.
    one_sided_phase : np.ndarray or None, optional
        One-sided phase array. If None, random phase will be generated, by default None.

    Returns
    -------
    np.ndarray
        Complex two-sided array.

    """
    if one_sided_phase is None:
        one_sided_phase = generate_one_sided_random_phase(one_sided_abs.size)
    double_sided_abs = one_sided_2_two_sided(one_sided_abs)
    double_sided_phase = one_sided_2_two_sided(one_sided_phase, sign=-1)
    return double_sided_abs * np.exp(1j * double_sided_phase)


def find_closest_index(arr: np.ndarray, value: float) -> int:
    """
    Find the index of the closest value in an array to a given value.

    Parameters
    ----------
    arr : np.ndarray
        Input array.
    value : float
        Target value.

    Returns
    -------
    int
        Index of the closest value.

    """
    return np.argmin(np.abs(arr - value))


def bin_df_rows(
    df: pd.DataFrame, binning_column: str, bins: Union[list, np.ndarray]
) -> pd.DataFrame:
    """
    Bin DataFrame rows based on a given column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    binning_column : str
        Column to bin on.
    bins : list or np.ndarray
        Bins to use for binning.

    Returns
    -------
    pd.DataFrame
        Binned DataFrame.

    """
    df["bins"] = pd.cut(df[binning_column], bins=bins, labels=bins[:-1])
    return df.groupby("bins").mean()


def abs_fft(
    timeTrace: Union[np.ndarray, pd.Series], fs: float = 250.0, window: float = 1.0
) -> np.ndarray:
    """
    Compute one-sided Fourier transform of a time series.

    Parameters
    ----------
    timeTrace : np.ndarray or pd.Series
        Time series to transform.
    fs : float, optional
        Sampling frequency, by default 250.0.
    window : float, optional
        Window to apply to the time series before transforming, by default 1.0.

    Returns
    -------
    np.ndarray
        One-sided Fourier transform.

    Raises
    ------
    ValueError
        If input time series is not a numpy array or pandas Series.

    """
    try:
        timeTrace = timeTrace.values
    except:
        pass
    spectrum = np.abs(fft(timeTrace))
    spectrum = spectrum[0 : int(timeTrace.size // 2) + 1]
    freq = np.arange(0, timeTrace.size // 2 + 1) * timeTrace.size // 2 * (fs / 2)
    return spectrum


class AppendOnlyIfMore(argparse.Action):
    """
    Custom action class for argparse that appends values to a list only if more than one value is provided.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The argument parser object.
    namespace : argparse.Namespace
        The namespace object where the parsed values are stored.
    values : List[str]
        The list of values provided for the option.
    option_string : str, optional
        The option string.

    Returns
    -------
    None

    Notes:
    ------
    - This action class is designed to be used with the argparse module.
    - It appends values to the list attribute specified by 'dest' only if more than one value is provided.
    - If only one value is provided, it replaces the list attribute with that single value.
    - It is intended for use with options that can accept multiple values.

    Example
    -------
    # Example usage
    parser = argparse.ArgumentParser()
    parser.add_argument("--values", action=AppendOnlyIfMore, default=[], dest="values", nargs="+")
    args = parser.parse_args()
    print(args.values)
    """

    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest, [])
        if items == self.default:
            items = []
        if len(values) == 1:
            items = values[0]
        else:
            items += values
        setattr(namespace, self.dest, items)


def distort_array(arr: np.ndarray, rstd: float = 0.5) -> np.ndarray:
    """
    Distort the values of a NumPy array by adding Gaussian noise.

    Parameters:
    -----------
    arr : np.ndarray
        The input NumPy array.
    rstd : float, optional
        The relative standard deviation of the Gaussian noise.

    Returns:
    --------
    np.ndarray
        The distorted array.

    Notes:
    ------
    - This function adds Gaussian noise to the values of the input array.
    - The amount of noise is controlled by the relative standard deviation (rstd).
    - A higher rstd value results in more distortion.
    - The function handles both masked and unmasked arrays.
    - For masked arrays, the distortion is only applied to unmasked values.

    Example:
    --------
    # Example usage
    arr = np.array([1, 2, 3, 4, 5])
    distorted_arr = distort_array(arr, rstd=0.5)
    print(distorted_arr)
    """
    if rstd == 0:
        return arr
    if isinstance(arr, np.ma.MaskedArray):
        mask = arr.mask
        arr = arr.data
        masked = True
    else:
        masked = False
    distorted_arr = np.asarray(
        [np.random.normal(mu, rstd * abs(mu), size=1)[0] for mu in arr]
    )

    if masked:
        distorted_arr = np.ma.array(distorted_arr, mask=mask)

    return distorted_arr


def autoscale_y(ax: axes.Axes, margin: float = 0.1) -> None:
    # modified: https://stackoverflow.com/questions/29461608/matplotlib-fixing-x-axis-scale-and-autoscale-y-axis
    """
    Rescales the y-axis based on the visible data given the current xlim of the axis. Credit: Dan Hickstein

    Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes object.
        margin (float): The fraction of the total height of the y-data to pad the upper and lower ylims.

    Returns:
        None
    """

    def get_bottom_top(line: np.ndarray) -> tuple[float, float]:
        """
        Calculates the bottom and top values of the y-data displayed within the x-limits.

        Parameters:
            line (np.ndarray): The y-data for a line.

        Returns:
            tuple[float, float]: The bottom and top values.
        """
        xd = line.get_xdata()
        yd = line.get_ydata()
        xd[np.abs(xd) == np.inf] = np.nan
        yd[np.abs(yd) == np.inf] = np.nan
        lo, hi = ax.get_xlim()
        y_displayed = yd[((xd > lo) & (xd < hi))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed) - margin * h
        top = np.max(y_displayed) + margin * h
        return bot, top

    lines = ax.get_lines()
    bot, top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot:
            bot = new_bot
        if new_top > top:
            top = new_top

    ax.set_ylim(bot, top)


def concatenate_simulated_dfs(dir_path: str, spec_substring: str = "") -> pd.DataFrame:
    """
    Concatenate simulated DataFrames from CSV files.

    Parameters
    ----------
    dir_path : str
        Directory path containing the CSV files.
    spec_substring : str, optional
        Substring to filter specific CSV files (default is an empty string).

    Returns
    -------
    pd.DataFrame
        Concatenated DataFrame containing the simulated data.

    """
    df_files = [
        os.path.join(dir_path, i)
        for i in os.listdir(dir_path)
        if (i.endswith(".csv") & (spec_substring in i))
    ]
    power_sim_DF = pd.read_csv(df_files[-1], index_col=0)
    # power_sim_DF.iloc[:, :] = power_sim_DF.values
    power_sim_DF.columns = power_sim_DF.columns.astype(float)

    df_list = []
    df_names = []
    # except_this = 'Salla_GSM16'
    except_this = "none"
    for f in df_files:
        if except_this not in f:
            df = pd.read_csv(f, index_col=0)
            df.columns = df.columns.astype(float)
            df_list.append(df)
            df_names.append(Path(f).stem)

    return pd.concat(df_list, keys=df_names)


def dropnans(arr: np.ndarray) -> np.ndarray:
    """
    Drop NaN values from the array.

    Parameters
    ----------
    arr : array-like
        Input array.

    Returns
    -------
    array
        Array with NaN values removed.

    """
    return arr[~np.isnan(arr)]
