#!/bin/bash -ue
trimmomatic PE -phred33     A_33_FDSW202661760-1r_HTNK5DSXY_L3_sub_1.fq.gz A_33_FDSW202661760-1r_HTNK5DSXY_L3_sub_2.fq.gz     A_33_FDSW202661760-1r_HTNK5DSXY_L3_sub_*.paired.fq.gz A_33_FDSW202661760-1r_HTNK5DSXY_L3_sub_*.unpaired.fq.gz     ILLUMINACLIP:adapters.fa:2:30:10
