#!/bin/bash
#-------
# USAGE
#-------
# sh PROTCLUS.sh                           # Running in foreground
# - OR -
# chmod 755 PROTCLUS.sh
# ./PROTCLUS.sh
# ------
# nohup sh PROTCLUS.sh > PROTCLUS.log &    # Running in background

#-------
# INTRO
#-------
# Created by Dheeraj Prakaash (2021)
#
# This script (PROTCLUS) was created to:
# 1. measure protein clustering vs time for a given trajectory (and plot it)
# 2. provide protein clustering score for a given trajectory
#
# Clustering score calculation: 
# If protein1 contacts 1 other protein, then the clustering score of protein1 = 1 for that frame
# If protein1 contacts 2 other proteins, then the clustering score of protein1 = 2 for that frame
# ... and so on
#
# Thus, PROTCLUS first calculates the clustering score of each protein per frame.
# Then, it calculates X = sum of clustering scores of each protein throughout the trajectory. 
# Finally, it adds X1 + X2 + X3 .... X'n' (n= num of proteins) to provide the final clustering score.
# Thus, the clustering score = sum of clustering scores of all proteins throughout the trajectory. 
#
# For multiple replicas, it is recommended to run this script for each trajectory (solvent and/or other unwanted components removed)

#--------------
# REQUIREMENTS
#--------------
# gromacs
# python3 
# numpy
# matplotlib

#-----------
# VARIABLES
#-----------
protein_count=16             # number of proteins in the system
protein_name="PROTEIN_NAME"
INDEX="group.ndx"            # index containing each protein as a separate group starting from group #1
TPR="../TPR_FILE.tpr"        # TPR of ONLY PROTEIN; Specify path if not in current dir
XTC="../XTC_FILE.xtc"        # XTC of ONLY PROTEIN; Specify path if not in current dir
TIME_BEGIN=0
TIME_END=24000000            # in picoseconds
TIME_PER_FRAME=60000         # in picoseconds
CONTACT_CUTOFF="0.55"        # in nanometers. Use "" for non-integers

#------------
# RUN SCRIPT
#------------

mkdir PROTCLUS_STEP_{1..3} 

# ------
# STEP 1
# ------

# For each pairwise combination of proteins
for (( i=1; i<=$protein_count; i++ ))
do
	for (( j=1; j<=$protein_count; j++ ))
	do
		if [ $i -ne $j ]
		then
			# run gmx mindist and move result files into corresponding directories
			mkdir PROTCLUS_STEP_1/$protein_name-${i}_with_${j}
			echo ${i} ${j} | gmx mindist -f $XTC -s $TPR -n $INDEX -od PROTCLUS_STEP_1/$protein_name-${i}_with_${j}/mindist_$protein_name-${i}_cont_with_${j}.xvg -or PROTCLUS_STEP_1/$protein_name-${i}_with_${j}/mindistres_$protein_name-${i}_cont_with_${j}.xvg -o PROTCLUS_STEP_1/$protein_name-${i}_with_${j}/atm-pair_$protein_name-${i}_cont_with_${j}.out -on PROTCLUS_STEP_1/$protein_name-${i}_with_${j}/numcont_$protein_name-${i}_cont_with_${j}.xvg -d $CONTACT_CUTOFF -pbc -respertime -printresname -b $TIME_BEGIN -e $TIME_END -group
		     
			# remove lines starting with # and @
			sed '/@/d;/#/d' PROTCLUS_STEP_1/$protein_name-${i}_with_${j}/numcont_$protein_name-${i}_cont_with_${j}.xvg > PROTCLUS_STEP_1/$protein_name-${i}_with_${j}/numcont_$protein_name-${i}_cont_with_${j}_clean.xvg
		     
			# extract column 2 of numcont...xvg
			awk '{print $2}' PROTCLUS_STEP_1/$protein_name-${i}_with_${j}/numcont_$protein_name-${i}_cont_with_${j}_clean.xvg > PROTCLUS_STEP_1/$protein_name-${i}_with_${j}/numcont_$protein_name-${i}_cont_with_${j}_col2.xvg

		fi
	done

# ------
# STEP 2 
# ------
	mkdir PROTCLUS_STEP_2/$protein_name-${i}
       
	# combine column 2 from all numcont...xvg files for each protein# 
	# ----------------
	# Example (prot#1):
	# -------------------------------------------------
	# (with prot#2)  (with prot#3)  (with prot#4) ...
	# value          value          value
	# .              .              .
	# .              .              .
	# -------------------------------------------------
	# do this for each prot#

	# Note: this does not paste columns in order (of prot#) unless $protein_count < 10 , but this does not affect calculating the sum of columns
	paste -d '\t' PROTCLUS_STEP_1/$protein_name-${i}_with_*/*col2.xvg > PROTCLUS_STEP_2/$protein_name-${i}/$protein_name-${i}_with_all_col2.xvg
	
	# replace non-zero values with 1
	awk 'BEGIN { OFS = FS = "\t" }
	{
	for (k = 1; k <= NF; ++k) 
	       { 
	       if ($k != "0") 
		       {
		       $k = "1";
		       }
	       }
	}
	{ print }' PROTCLUS_STEP_2/$protein_name-${i}/$protein_name-${i}_with_all_col2.xvg > PROTCLUS_STEP_2/$protein_name-${i}/$protein_name-${i}_with_all_col2_nz1.xvg
	
	# generate sum of all contacts (columns) of a prot# with all other prot# per timestep (row)
	awk '
	{
	sum = 0; 
	for (m = 1; m <= NF; m++) 
		{
		sum += $m; 
		}
	print sum
	}' PROTCLUS_STEP_2/$protein_name-${i}/$protein_name-${i}_with_all_col2_nz1.xvg > PROTCLUS_STEP_2/$protein_name-${i}/$protein_name-${i}_col2_sum.xvg
	
	# transpose column to row for each prot#
	paste -s PROTCLUS_STEP_2/$protein_name-${i}/$protein_name-${i}_col2_sum.xvg > PROTCLUS_STEP_2/$protein_name-${i}/$protein_name-${i}_col2_sum_tr.xvg
done

# ------
# STEP 3
# ------
# (1) calc total contacts per timestep
# (2) calc total contacts per protein
# (3) calc CLUST_SCORE = Sum of total contacts per prot
# (4) prepare final dataframe 
# --------------------------------------------------
# 0  1  2  3 ... x $TIME_PER_FRAME (ps)
# 1  .  .  .
# 2  .  .  .
# 3  .  .  .
# |  .  .  .
# v  .  .  .
# prot#
# --------------------------------------------------

eval cat PROTCLUS_STEP_2/$protein_name-{1..$protein_count}/*tr.xvg > PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-raw_df.xvg

# calc sum of rows (i.e. per frame) and columns (i.e. per protein) 

awk '
{ 
TR=0
for ( n = 1; n <= NF; n++ ) 
{
	TR += $n
	TC[n] += $n
	printf( "%6d", $n )
}
print "  = " TR
TF = NF
}
END {
for ( n = 1; n <= NF; n++ )
	printf "%6d", TC[n]
	print ""
}' PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-raw_df.xvg > PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-sum_of_cols_and_rows.xvg


#---(1) TOTAL CONTACTS PER TIMESTEP ---

tail -n 1 PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-sum_of_cols_and_rows.xvg > PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-cont_sum_per_timestep.xvg

tr -s "[:space:]" "\n" < PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-cont_sum_per_timestep.xvg | tail -n +2 > PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-cont_sum_per_timestep_tr.xvg

seq $TIME_BEGIN $TIME_PER_FRAME $TIME_END > PROTCLUS_STEP_2/timesteps_col.xvg

paste -d '\t' PROTCLUS_STEP_2/timesteps_col.xvg PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-cont_sum_per_timestep_tr.xvg > PROTCLUS_STEP_3/$protein_name-1_to_$protein_count-TOTAL_CONT_PER_TIMESTEP.xvg


#---(2) TOTAL CONTACTS PER PROTEIN ---

awk '{print $NF}' PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-sum_of_cols_and_rows.xvg | head -n -1 > PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-cont_sum_per_prot.xvg

seq 1 $protein_count > PROTCLUS_STEP_2/protnum_col_1-$protein_count.xvg

paste -d '\t' PROTCLUS_STEP_2/protnum_col_1-$protein_count.xvg PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-cont_sum_per_prot.xvg > PROTCLUS_STEP_3/$protein_name-1_to_$protein_count-TOTAL_CONT_PER_PROTEIN.xvg


#---(3) CLUSTERING SCORE ---

paste -s PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-cont_sum_per_prot.xvg > PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-cont_sum_per_prot_tr.xvg

awk '
{
sum = 0; 
for (p = 1; p <= NF; p++) 
	{
	sum += $p; 
	}
print sum
}' PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-cont_sum_per_prot_tr.xvg > PROTCLUS_STEP_3/$protein_name-1_to_$protein_count-CLUST_SCORE.xvg


#---(4) PROTCLUS vs TIME DATAFRAME ---

seq -s "	" $TIME_BEGIN $TIME_PER_FRAME $TIME_END > PROTCLUS_STEP_2/timesteps_row.xvg

seq 0 $protein_count > PROTCLUS_STEP_2/protnum_col_0-$protein_count.xvg

eval cat PROTCLUS_STEP_2/timesteps_row.xvg PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-raw_df.xvg > PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-df_with_timesteps.xvg

paste -d '\t' PROTCLUS_STEP_2/protnum_col_0-$protein_count.xvg PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-df_with_timesteps.xvg > PROTCLUS_STEP_3/$protein_name-1_to_$protein_count-PROTCLUS_vs_TIME_DATAFRAME.xvg


# --------
# PLOTTING
# --------
cp PROTCLUS_STEP_2/$protein_name-1_to_$protein_count-raw_df.xvg PROTCLUS_STEP_3/raw_data.xvg
cd PROTCLUS_STEP_3

echo "import matplotlib.pyplot as plt" >> plot_PROTCLUS_vs_TIME.py
echo "import numpy as np" >> plot_PROTCLUS_vs_TIME.py
echo "data = np.loadtxt('raw_data.xvg')" >> plot_PROTCLUS_vs_TIME.py
echo "fig, ax = plt.subplots()" >> plot_PROTCLUS_vs_TIME.py
echo "im = ax.pcolormesh(data, cmap='RdPu')" >> plot_PROTCLUS_vs_TIME.py
echo "plt.colorbar(im, label='Number of proteins in contact')" >> plot_PROTCLUS_vs_TIME.py
echo "plt.xlabel('Frames', fontsize=13)" >> plot_PROTCLUS_vs_TIME.py
echo "plt.ylabel('Protein index', fontsize=13)" >> plot_PROTCLUS_vs_TIME.py
echo "fig.savefig('raw_data.svg', format='svg', dpi=600)" >> plot_PROTCLUS_vs_TIME.py
echo "fig.savefig('raw_data.png', dpi=1200)" >> plot_PROTCLUS_vs_TIME.py
python3 plot_PROTCLUS_vs_TIME.py
