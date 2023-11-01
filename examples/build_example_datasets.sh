#!/bin/bash

GMAPS='LFSS GSM08 GSM16 Haslam LFmap SSM GMOSS ULSA'

#base=./
base=./examples

: <<'END'
python ${base}/utilities/create_isotropic_antenna_pattern.py -o ${base}/antenna_setup_files/Isotropic_antenna_pattern.xml
python ${base}/utilities/create_flat_hw_profile.py -o ${base}/antenna_setup_files/HardwareProfileList_flat.xml

###  
echo **************** simulated power in HW ****************
# realistic
for GMAP in $GMAPS; do
    echo $GMAP;
    python ${base}/utilities/get_power_dataset.py -g $GMAP -a ${base}/antenna_setup_files/SALLA_EW.xml -w ${base}/antenna_setup_files/HardwareProfileList_realistic.xml -s Salla_EW -o ${base}/simulated_power_datasets/ -r EW
    echo =================================================================
done

#for GMAP in $GMAPS; do
#    echo $GMAP;
#    python ${base}/utilities/get_power_dataset.py -g $GMAP -a ${base}/antenna_setup_files/SALLA_NS.xml -w ${base}/#antenna_setup_files/HardwareProfileList_realistic.xml -s Salla_NS -o ${base}/simulated_power_datasets/ -r NS -z -180
#    echo =================================================================
#done

# idealistic
for GMAP in $GMAPS; do
    echo $GMAP;
    python ${base}/utilities/get_power_dataset.py -g $GMAP -a ${base}/antenna_setup_files/Isotropic_antenna_pattern.xml -w ${base}/antenna_setup_files/HardwareProfileList_flat.xml -s isoAnt_flathw -o ${base}/simulated_power_datasets/;
    echo =================================================================
done


### simulated voltage square density
echo **************** simulated voltage square density ****************
# realistic
for GMAP in $GMAPS; do
    echo $GMAP;
    python ${base}/utilities/get_power_dataset.py -g $GMAP -a ${base}/antenna_setup_files/SALLA_EW.xml -w -s Salla_EW -o ${base}/voltage2_density/ -t voltage2 -f 10 125;
    echo =================================================================
done

#idealistic    
for GMAP in $GMAPS; do
    echo $GMAP;
    python ${base}/utilities/get_power_dataset.py -g $GMAP -a ${base}/antenna_setup_files/Isotropic_antenna_pattern.xml -s isoAnt -o ${base}/voltage2_density/ -t voltage2 -f 10 125;
    echo =================================================================
done

# mock power datasets
echo **************** mock power datasets ****************
for GMAP in $GMAPS; do
    echo $GMAP;
python ${base}/utilities/get_mock_dataset.py -v ${base}/voltage2_density/voltage2_density_Salla_EW_$GMAP.csv -w ${base}/antenna_setup_files/HardwareProfileList_realistic.xml -n 10000 -a 0 -o ${base}/mock_power_datasets/
    echo =================================================================
done

for GMAP in $GMAPS; do
    echo $GMAP;
python ${base}/utilities/get_mock_dataset.py -v ${base}/voltage2_density/voltage2_density_isoAnt_$GMAP.csv -w ${base}/antenna_setup_files/HardwareProfileList_flat.xml -n 10000 -a 0 -o ${base}/mock_power_datasets/
    echo =================================================================
done

# various size
python ${base}/utilities/get_mock_dataset.py -v ${base}/voltage2_density/voltage2_density_Salla_EW_GSM16.csv -w ${base}/antenna_setup_files/HardwareProfileList_realistic.xml -n 5000 -a 0 -tor -o ${base}/mock_power_datasets/
python ${base}/utilities/get_mock_dataset.py -v ${base}/voltage2_density/voltage2_density_Salla_EW_GSM16.csv -w ${base}/antenna_setup_files/HardwareProfileList_realistic.xml -n 1000 -a 0 -tor -o ${base}/mock_power_datasets/

python ${base}/utilities/get_mock_dataset.py -v ${base}/voltage2_density/voltage2_density_Salla_EW_GSM16.csv -w ${base}/antenna_setup_files/HardwareProfileList_realistic.xml -n 5000 -a 0 -o ${base}/mock_power_datasets/
python ${base}/utilities/get_mock_dataset.py -v ${base}/voltage2_density/voltage2_density_Salla_EW_GSM16.csv -w ${base}/antenna_setup_files/HardwareProfileList_realistic.xml -n 1000 -a 0 -o ${base}/mock_power_datasets/

### Constant temperature ###
## without rounding, 
python ${base}/utilities/get_mock_dataset.py -v ${base}/voltage2_density/voltage2_density_Salla_EW_GSM16.csv -w ${base}/antenna_setup_files/HardwareProfileList_realistic.xml -n 10000 -a 0 -tor -t 30 -o ${base}/mock_power_datasets/
python ${base}/utilities/get_mock_dataset.py -v ${base}/voltage2_density/voltage2_density_isoAnt_GSM16.csv -w ${base}/antenna_setup_files/HardwareProfileList_flat.xml -n 10000 -a 0 -tor -t 30 -o ${base}/mock_power_datasets/

# with extra noise
python ${base}/utilities/get_mock_dataset.py -v ${base}/voltage2_density/voltage2_density_isoAnt_GSM16.csv -w ${base}/antenna_setup_files/HardwareProfileList_flat.xml -n 10000 -a 5e-16 -tor -t 30 -o ${base}/mock_power_datasets/


## with rounding 
python ${base}/utilities/get_mock_dataset.py -v ${base}/voltage2_density/voltage2_density_Salla_EW_GSM16.csv -w ${base}/antenna_setup_files/HardwareProfileList_realistic.xml -n 10000 -a 0 -t 30 -o ${base}/mock_power_datasets/

python ${base}/utilities/get_mock_dataset.py -v ${base}/voltage2_density/voltage2_density_isoAnt_GSM16.csv -w ${base}/antenna_setup_files/HardwareProfileList_flat.xml -n 10000 -a 0 -t 30 -o ${base}/mock_power_datasets/
END


# python convertSfile2xml.py -i ${base}/Sparamaters/RD_mean_LNA_EW.s2p -s HardwareProfileList_realistic.xml -m w -l LNA
# python convertSfile2xml.py -i ${base}/Sparamaters/Cable_assembly_ant-dig_with_pigtail.s2p -s HardwareProfileList_realistic.xml -m a -l cable_fromLNA2digitizer
# python convertSfile2xml.py -i ${base}/Sparamaters/Digitizer_v3_high_gain_MagAngle_dB_deg.s2p -s HardwareProfileList_realistic.xml -m a -l digitizer





