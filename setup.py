#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from setuptools.command.install import install
import urllib.request
import tarfile
import os
import sys
import shutil

# import wget


def readme():
    with open("README.md") as f:
        return f.read()


# prepare progressbar
def show_progress(block_num, block_size, total_size):
    print(round(block_num * block_size / total_size * 100, 2), end="\r")


# Define the url of the tar.gz package and the name of the folder inside it
package_url = "http://tianlai.bao.ac.cn/~huangqizhi/SSM.tar.gz"
package_file_name = "SSM.tar.gz"
package_folder = "SSM"

# Define a custom installation command that will download and extract the package before installing it
class ssm_install(install):
    def run(self):
        if not os.path.exists(package_file_name):
            try:
                # Import tqdm if it's installed
                from tqdm import tqdm

                # Download the package and show a progress bar
                with tqdm(
                    unit="B", unit_scale=True, desc=package_url.split("/")[-1]
                ) as progress:
                    package_file, _ = urllib.request.urlretrieve(
                        package_url,
                        reporthook=lambda blocknum, blocksize, totalsize: progress.update(
                            blocknum * blocksize - progress.n
                        ),
                    )
            except ImportError:
                # tqdm is not installed, download the package without a progress bar
                package_file, _ = urllib.request.urlretrieve(package_url)
        else:
            # The package has already been downloaded, use the existing file
            package_file = package_file_name
        # Extract the package to a temporary directory
        print(package_file)
        with tarfile.open(package_file) as tar:
            tar.extractall()
        # Call the original install command
        install.do_egg_install(self)
        shutil.rmtree(package_folder)


required = [
    "numpy",
    "scipy",
    "healpy",
    "pandas",
    "lxml",
    "statsmodels",
    "tqdm",
    "plotly",
    "pygdsm",
    "wget",
    "pygdsm",
    "pylfmap",
    "sklearn"
]


setup(
    name="radiocalibrationtoolkit",
    version="0.1",
    description="Tools for the absolute antenna calibration using the galactic radio emission",
    long_description=readme(),
    author="T. Fodran",
    author_email="t.fodran@astro.ru.nl",
    url="https://github.com/F-Tomas/radiocalibrationtoolkit.git",
    packages=find_packages(),
    install_requires=required,
    setup_requires=required,
    cmdclass={
        'install': ssm_install,
    },
    extras_require={
        "jupyter": ["jupyter"],
    },
    package_data={
        "radiocalibrationtoolkit": [
            "mock_power_datasets/*",
            "skymapwrappers/tabulated/*",
        ]
    },
)
