#!/bin/bash

GMAPS='LFSS GSM08 GSM16 Haslam LFmap SSM GMOSS ULSA'

python ./utilities/create_isotropic_antenna_pattern.py
python ./utilities/create_flat_hw_profile.py

# simulated power in HW
for GMAP in $GMAPS; do
    echo $GMAP;
    python ./utilities/get_power_dataset.py -g $GMAP -a ./antenna_setup_files/SALLA_EW.xml -w ./antenna_setup_files/HardwareProfileList_realistic.xml -s Salla;
    echo =================================================================
done
    
for GMAP in $GMAPS; do
    echo $GMAP;
    python ./utilities/get_power_dataset.py -g $GMAP -a ./antenna_setup_files/Isotropic_antenna_pattern.xml -w ./antenna_setup_files/HardwareProfileList_flat.xml -s isoAnt_flathw;
    echo =================================================================
done

# simulated voltage square density
for GMAP in $GMAPS; do
    echo $GMAP;
    python ./utilities/get_power_dataset.py -g $GMAP -a ./antenna_setup_files/SALLA_EW.xml -w -s Salla_EW -o ./voltage2_density/ -t voltage2 -f 10 125;
    echo =================================================================
done
    
for GMAP in $GMAPS; do
    echo $GMAP;
    python ./utilities/get_power_dataset.py -g $GMAP -a ./antenna_setup_files/Isotropic_antenna_pattern.xml -s isoAnt -o ./voltage2_density/ -t voltage2 -f 10 125;
    echo =================================================================
done

# mock power datasets
for GMAP in $GMAPS; do
    echo $GMAP;
python ./utilities/get_mock_dataset.py -v ./voltage2_density/voltage2_density_Salla_EW_$GMAP.csv -w ./antenna_setup_files/HardwareProfileList_realistic.xml -n 10000 -a 0
    echo =================================================================
done

for GMAP in $GMAPS; do
    echo $GMAP;
python ./utilities/get_mock_dataset.py -v ./voltage2_density/voltage2_density_isoAnt_$GMAP.csv -w ./antenna_setup_files/HardwareProfileList_flat.xml -n 10000 -a 0
    echo =================================================================
done

# various size
python ./utilities/get_mock_dataset.py -v ./voltage2_density/voltage2_density_Salla_EW_GSM16.csv -w ./antenna_setup_files/HardwareProfileList_realistic.xml -n 5000 -a 0 -tor
python ./utilities/get_mock_dataset.py -v ./voltage2_density/voltage2_density_Salla_EW_GSM16.csv -w ./antenna_setup_files/HardwareProfileList_realistic.xml -n 1000 -a 0 -tor

python ./utilities/get_mock_dataset.py -v ./voltage2_density/voltage2_density_Salla_EW_GSM16.csv -w ./antenna_setup_files/HardwareProfileList_realistic.xml -n 5000 -a 0
python ./utilities/get_mock_dataset.py -v ./voltage2_density/voltage2_density_Salla_EW_GSM16.csv -w ./antenna_setup_files/HardwareProfileList_realistic.xml -n 1000 -a 0

### Constant temperature ###
## without rounding, 
python ./utilities/get_mock_dataset.py -v ./voltage2_density/voltage2_density_Salla_EW_GSM16.csv -w ./antenna_setup_files/HardwareProfileList_realistic.xml -n 10000 -a 0 -tor -t 30
python ./utilities/get_mock_dataset.py -v ./voltage2_density/voltage2_density_isoAnt_GSM16.csv -w ./antenna_setup_files/HardwareProfileList_flat.xml -n 10000 -a 0 -tor -t 30

# with extra noise
python ./utilities/get_mock_dataset.py -v ./voltage2_density/voltage2_density_isoAnt_GSM16.csv -w ./antenna_setup_files/HardwareProfileList_flat.xml -n 10000 -a 5e-16 -tor -t 30


## with rounding 
python ./utilities/get_mock_dataset.py -v ./voltage2_density/voltage2_density_Salla_EW_GSM16.csv -w ./antenna_setup_files/HardwareProfileList_realistic.xml -n 10000 -a 0 -t 30
python ./utilities/get_mock_dataset.py -v ./voltage2_density/voltage2_density_isoAnt_GSM16.csv -w ./antenna_setup_files/HardwareProfileList_flat.xml -n 10000 -a 0 -t 30



# python convertSfile2xml.py -i ./Sparamaters/RD_mean_LNA_EW.s2p -s HardwareProfileList_realistic.xml -m w -l LNA
# python convertSfile2xml.py -i ./Sparamaters/Cable_assembly_ant-dig_with_pigtail.s2p -s HardwareProfileList_realistic.xml -m a -l cable_fromLNA2digitizer
# python convertSfile2xml.py -i ./Sparamaters/Digitizer_v3_high_gain_MagAngle_dB_deg.s2p -s HardwareProfileList_realistic.xml -m a -l digitizer





