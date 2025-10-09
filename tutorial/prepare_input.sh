#!/bin/bash

# input filename
pr0max_dark_model="./input/8Z1J.pdb"
Fodark_file="./input/8Z1J.mtz"
Fotr_file="./input/8Z3X.mtz"

# input file parameters
low_resolution=30.00
high_resolution=1.64
symmetry="'P 212121'"

# Check if the three input file exists
if [[ ! -f $pr0max_dark_model ]] ; then
    echo "FileNotFoundError: File $pr0max_dark_model not found."
    exit
fi
if [[ ! -f $Fodark_file ]] ; then
    echo "FileNotFoundError: File $Fodark_file not found."
    exit
fi
if [[ ! -f $Fotr_file ]] ; then
    echo "FileNotFoundError: File $Fotr_file not found."
    exit
fi

echo "Start preparing the input file for dFoCC..."
echo "If the script fails in the middle, please check the log file: ./input/prepare_input.log"

{
    # Check if any command returns a non-zero value (fails)
    set -e

    FoFcdark_file="${Fodark_file%.mtz}-fofc.mtz"

    sfall XYZIN $pr0max_dark_model HKLIN $Fodark_file HKLOUT $FoFcdark_file <<EOF
TITLE Calculate Dark Structure Factor
LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag
LABOUT FC=Fcdark PHIC=PHIcdark
MODE SFCALC XYZIN HKLIN
RESOLUTION $low_resolution $high_resolution
SYMMETRY $symmetry
END
EOF

    # Combine dark and triggered structure factor and rename several columns
    FoFcdark_and_Fotr_file="${Fotr_file%.mtz}"-combined.mtz

    cad HKLIN1 $FoFcdark_file HKLIN2 $Fotr_file HKLOUT $FoFcdark_and_Fotr_file <<EOF
TITLE Combine Dark and Triggered Structure Factor
MONITOR BRIEF
LABIN file 1 E1=FreeR_flag E2=FP E3=SIGFP E4=Fcdark E5=PHIcdark
LABOUT file 1 E1=FreeR_flag E2=Fodark E3=SIGFodark E4=Fcdark E5=PHIcdark
CTYPIN file 1 E1=I E2=F E3=Q E4=F E5=P
RESOLUTION file 1 $low_resolution $high_resolution
LABIN file 2 E1=FP E2=SIGFP
LABOUT file 2 E1=Fotr E2=SIGFotr
CTYPIN file 2 E1=F E2=Q
RESOLUTION file 2 $low_resolution $high_resolution
EOF

#     # Alternatively, combine them with `sftools`
#     sftools <<EOF
# read $FoFcdark_file col FreeR_flag FP SIGFP Fcdark PHIcdark
# set label col FP
# Fodark
# set label col SIGFP
# SIGFodark
# read $Fotr_file col FP SIGFP
# set label col FP
# Fotr
# set label col SIGFP
# SIGFotr
# write $FoFcdark_and_Fotr_file
# end
# EOF

    # Scale the Fotr dataset against Fodark
    FoFcdark_and_Fotr_scaled_file="${Fotr_file%.mtz}-scaled.mtz"

    scaleit HKLIN $FoFcdark_and_Fotr_file HKLOUT $FoFcdark_and_Fotr_scaled_file <<EOF
TITLE Scale Triggered Structure Factor
RESOLUTION 10.0 2.5
NOWT
CONVERGE NCYC 4
CONVERGE ABS 0.001
CONVERGE TOLR -7
REFINE SCALE WILSON
AUTO
LABIN FP=Fodark SIGFP=SIGFodark FPH1=Fotr SIGFPH1=SIGFotr
EOF

    # Split the dark and the triggered dataset, so that Phenix could read them correctly
    Fodark_split_file="${Fodark_file%.mtz}-split.mtz"
    Fotr_split_file="${Fotr_file%.mtz}-split.mtz"

    cad HKLIN1 $FoFcdark_and_Fotr_scaled_file HKLOUT $Fodark_split_file <<EOF
TITLE Split Fodark
MONITOR BRIEF
LABIN file 1 E1=FreeR_flag E2=Fodark E3=SIGFodark
LABOUT file 1 E1=FreeR_flag E2=Fodark E3=SIGFodark
CTYPIN file 1 E1=I E2=F E3=Q
RESOLUTION file 1 $low_resolution $high_resolution
EOF

    cad HKLIN1 $FoFcdark_and_Fotr_scaled_file HKLOUT $Fotr_split_file <<EOF
TITLE Split Fotr
MONITOR BRIEF
LABIN file 1 E1=Fotr E2=SIGFotr
LABOUT file 1 E1=Fotr E2=SIGFotr
CTYPIN file 1 E1=F E2=Q
RESOLUTION file 1 $low_resolution $high_resolution
EOF

# # Alternatively, split them with `sftools`
# sftools <<EOF
# read $FoFcdark_and_Fotr_scaled_file
# write $Fodark_split_file col FreeR_flag Fodark Fcdark
# write $Fotr_split_file col Fotr SIGFotr
# END
# EOF

    # Generate dFo structure factor
    dFo_file="${Fotr_file%.mtz}-dFo.mtz"

    phenix.fobs_minus_fobs_map job_title="Generate dFo Structure Factor" \
        f_obs_1_file=$Fotr_split_file \
        f_obs_2_file=$Fodark_split_file \
        output_file=$dFo_file \
        f_obs_1_label=Fotr \
        f_obs_2_label=Fodark \
        phase_source=$pr0max_dark_model \
        high_res=$high_resolution low_res=$low_resolution \
        sigma_cutoff=3.0
} > ./input/prepare_input.log

echo "File preparation is done."
echo "Please modify the configuration file (./configuration.phil) with the prepared file:"
echo "input.pdb.dark_file_name = ${pr0max_dark_model}"
echo "input.dark_data.filename = ${FoFcdark_and_Fotr_scaled_file}"
echo "input.difference_data.filename = ${dFo_file}"
