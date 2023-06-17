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

The ULSA is currently tabulated only in the 1-125 MHz range for frequency dependant index. The GMOSS up to 400 MHz.

GMOSS, Haslam, and GSM2008 are without CMB. To fix this issue, a value of 2.7255 is added to GMOSS upon reading the tabulated values. The Haslam and GSM2008 are fixed for the CMB using a decorator around their `pygdsm` classes. So, do not import them manually from `pygdsm`; their fixed versions are imported by 

`from radiocalibrationtoolkit import *`

All the example's data files are pre-generated, so each example should work out of the box. Also, all of these data files can be reproduced by running the scripts.


## Install

run this in terminal:

`python -m pip install git+https://github.com/F-Tomas/radiocalibrationtoolkit.git`

## Documentation

Check the documentation and tutorial [here](https://f-tomas.github.io/radiocalibrationtoolkit/index.html)

To compile the documentation manually, go to `docs` and execute: 

`sphinx-build -b html . ./html`

