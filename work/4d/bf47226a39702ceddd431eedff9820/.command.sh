#!/bin/bash -ue
trimmomatic PE -phred33 A_33_FDSW202661760-1r_HTNK5DSXY_L3_sub_1.fq.gz A_33_FDSW202661760-1r_HTNK5DSXY_L3_sub_2.fq.gz A_33_FDSW202661760-1r_HTNK5DSXY_L3_sub_1.trimmed.fq.gz A_33_FDSW202661760-1r_HTNK5DSXY_L3_sub_1.discarded.fq.gz A_33_FDSW202661760-1r_HTNK5DSXY_L3_sub_2.trimmed.fq.gz A_33_FDSW202661760-1r_HTNK5DSXY_L3_sub_2.discarded.fq.gz 
ILLUMINACLIP:/Users/vinithanadar/miniconda/envs/nf_course/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10
