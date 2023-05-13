# Radio calibration toolkit


## Getting started

Temporary version!

The package uses 5 different galactic model softwares:

[SSM](http://tianlai.bao.ac.cn/~huangqizhi/)

[GMOSS](https://github.com/mayurisrao/GMOSS)

[ULSA](https://github.com/Yanping-Cong/ULSA/tree/v2.0)

[GSM2016, GSM2008, LFSS, Haslam](https://github.com/telegraphic/pygdsm)

[LFmap](http://www.astro.umd.edu/~emilp/LFmap/LFmap_1.0.tar) using [pylfmap](https://github.com/F-Tomas/pylfmap)


The GMOSS maps are tabulated, because the code is in C++, ULSA takes very long to generate one map, so I pre-generate them and tabulated them. The SSM and pyGDSM are fast and in python, so no need to tabulate these maps.

The LFmap is in C++, but I wrote a wrapper/interface so that the model can be imported as a Python package.

Thus, you need to install pygdsm (pip install pygdsm), the SSM and ULSA can be installed by cloning them and using the setup.py, the pylfmap is ready to go, you just need to add it to your python path enviromental variable.

The ULSA is currently tabulated only in 30-80MHz for frequency dependant index. Likewise the GMOSS. I will generate tables for these for the full 0-125 range (if possible).

All the examples data files are pre-generated, so each example should work out of the box. Also, all of these data files can be reproduce by running the scripts.


## Install

run this in terminal:

`python -m pip install git+https://github.com/F-Tomas/radiocalibrationtoolkit.git`

## Documentation

Check the documentation and tutorial [here](https://f-tomas.github.io/radiocalibrationtoolkit/index.html)

To compile the documentation manually, go to `docs` and execute: 

`sphinx-build -b html . ./html`

