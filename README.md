# Radio calibration toolkit


## Getting started

The package uses 5 different galactic model softwares:

[SSM](http://tianlai.bao.ac.cn/~huangqizhi/)

[GMOSS](https://github.com/mayurisrao/GMOSS)

[ULSA](https://github.com/Yanping-Cong/ULSA/tree/v2.0)

[GSM2016, GSM2008, LFSS, Haslam](https://github.com/telegraphic/pygdsm)

[LFmap](http://www.astro.umd.edu/~emilp/LFmap/LFmap_1.0.tar) using [pylfmap](https://github.com/F-Tomas/pylfmap)


The GMOSS maps are tabulated because the code is in C++, ULSA takes very long to generate one map, so I pre-generate and tabulate them. The SSM and pyGDSM are fast and in Python, so tabulating these maps is unnecessary.

The LFmap is in C++, but I wrote a wrapper/interface to import the model as a Python package.

Thus, you need to install pygdsm (pip install pygdsm); the SSM and ULSA can be installed by cloning them and using the setup.py. The pylfmap is now available at pip: `pip install pylfmap`.

The ULSA maps with frequency dependent index and GMOSS maps are tabulated up to 400 MHz in 1 MHz steps.

GMOSS, Haslam, and GSM2008 are without CMB. To fix this issue, a value of 2.7255 is added to GMOSS upon reading the tabulated values. The Haslam and GSM2008 are fixed for the CMB using a decorator around their `pygdsm` classes. So, do not import them manually from `pygdsm`; their fixed versions are imported by 

`from radiocalibrationtoolkit import *`

All the example's data files are pre-generated, so each example should work out of the box. Also, all of these data files can be reproduced by running the scripts.


## Install

run this in terminal:

`python -m pip install git+https://github.com/F-Tomas/radiocalibrationtoolkit.git`

## Documentation

Documentation and tutorial for version `v0.1.0-alpha` is available [here](https://f-tomas.github.io/radiocalibrationtoolkit/index.html)

To compile the documentation manually, go to `docs` and execute: 

`sphinx-build -b html . ./html`


## Some tips on installing the ULSA

Installing the [ULSA](https://github.com/Yanping-Cong/ULSA/tree/v2.0) package was not quite straightforward for me. Hence I share briefly some of my experience from the installation.

First:
- download and install: l_fortran-compiler_p_XXXX.X.X.XXXXX_offline.sh from [here]( https://registrationcenter-download.intel.com/akdlm/IRC_NAS/150e0430-63df-48a0-8469-ecebff0a1858/). I worked with this version: l_fortran-compiler_p_2023.0.0.25394_offline.sh. Once installed, `source setvars.sh` from the `oneapi` directory,

- clone [caput](https://github.com/zuoshifan/caput.git) and switch to branch `origin/zuo/develop` and after fixing it to work with `Python 3` install it,
 
- to make ULSA work with `Python 3`, it is good to install [2to3](https://docs.python.org/3/library/2to3.html) convertor and run it on all `.py` files in ULSA and caput. Some changes you will need to do manually, for example, exchange `np.int` for just `int` and fix slicing of arrays in file `mpiutil.py` from `[start:stop]` to `[int(start):int(stop)]`,

- install all other required dependencies.

To create the `libNE2001.so` file, once unpacked NE2001_4python.zip go to src.NE2001 directory and in all `.f` files, replace the relative path with the absolute path to the required files.

E.g. in file `nevoidN.f`, change 
```
	  open(luvoid, file='nevoidN.NE2001.dat', status='old')
```
to
```
	  open(luvoid, file='/vol/astro6/auger-radiodigitizer/'
     &//'skymaps/ULSA/ULSA/NE2001/NE2001_4python/'
     &//'NE2001_4python/bin_NE2001/nevoidN.NE2001.dat',
     &status='old')
```
and so on. The `&//` is the Fortran way of splitting the path into multiple lines. Once fixed, you can run `make so` to generate the `libNE2001.so` file. Then find all references to this file in ULSA `.py` files and rewrite them to your path.

Install ULSA, and it should work now.

