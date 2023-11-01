import os
import sys

from pylfmap import LFmap

from radiocalibrationtoolkit import *

LATITUDE = -35.206667

antenna_inst = AntennaPattern(
    "./antenna_setup_files/SALLA_32S_4NEC2_WCD_8S_PERMITTIVITY5.5_CONDUCTIVITY0.0014_extended_EW.xml",
)

# Prepare objects
lfmap = LFmap()
lfss = LowFrequencySkyModel(freq_unit="MHz")
gsm2008 = GlobalSkyModel(freq_unit="MHz")
gsm2016 = GlobalSkyModel2016(freq_unit="MHz")
haslam = HaslamSkyModel(freq_unit="MHz", spectral_index=-2.53)
ssm = SSM()
gmoss = GMOSS()
ulsa_fdi = ULSA(index_type="freq_dependent_index")


map_instances_dict = {
    "LFSS": lfss,
    "GSM08": gsm2008,
    "GSM16": gsm2016,
    "Haslam": haslam,
    "LFmap": lfmap,
    "SSM": ssm,
    "GMOSS": gmoss,
    "ULSA": ulsa_fdi,
}

# for just the voltage square spectral density
# do not use impedance and do not integrate!

for map_label in tqdm(map_instances_dict.keys()):
    galactic_map_inst = map_instances_dict[map_label]

    lst_range = range(25)
    freq_Mhz_range = range(10, 125, 1)

    voltage2_density_DF = calculate_power_spectral_density(
        antenna_inst=antenna_inst,
        galactic_map_inst=galactic_map_inst,
        lst_range=lst_range,
        freq_Mhz_range=freq_Mhz_range,
        latitude=LATITUDE,
        update_antenna_conventions={
            "shift_phi": -90,
            "flip_theta": True,
            "flip_phi": False,
            "in_degrees": True,
            "add_invisible_sky": True,
        },
    )
    voltage2_density_DF.index.name = "LST"
    voltage2_density_DF.columns.name = "freq_MHz"

    voltage2_density_DF.to_csv(
        "./voltage2_density/voltage2_density_{}.csv".format(map_label)
    )
