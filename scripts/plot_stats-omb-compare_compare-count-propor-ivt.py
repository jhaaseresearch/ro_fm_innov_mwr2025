#________________________________________________________________
#
#   Python script to read in ERA5 (Renalyses) in
#    CIRAQ-style and ASCII format and
#    occultation atmPrf netCDF files and plot the
#    difference between them.
#   This script calculates the stats for all occs
#    from all IOPs for all years for the given
#    dataset.
#
#   ~ mjm 2020-09-02
#________________________________________________________________

import numpy as np
import matplotlib
matplotlib.use('pdf')  # required for interactive backend node
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import subprocess
import sys
import os
from datetime import datetime, timedelta
sys.path.insert(0, './calc_vars')
from label_axes_iterate import label_axes

# Read args
nargs = len(sys.argv) - 1
script_name = sys.argv[0]  # this will always work.

nargs_req = 4
if nargs == nargs_req:
	arrecon_year_str = sys.argv[1]
	iop_str = sys.argv[2]
	series_id = sys.argv[3]
	io_root = sys.argv[4]
else:
	print("")
	print("")
	print("   Wrong number of arguments, ABORTING!!")
	print("")
	print("   You gave " + str(nargs) + " of the " + str(nargs_req) + " required args")
	print("")
	print("     Usage:")
	print("       " + script_name + " <arrecon_year> <iop_str> <series_id> <output_root>")
	print("")
	print("     Example:")
	print("       python " + script_name + " 2022 full 05 ../")
	print("")
	print("")

	sys.exit()

# User input
inc_hgt_bins_m = 200

missingValue = -999.0  # this same missing value seems to be used for all variables in soundings

nlines_header = 1

input_dir_root = io_root
output_dir_root = io_root
output_path_pdf = output_dir_root + 'plots/'

# Make list of input files
PATH_OBS = input_dir_root + 'stats_ascii/'
FILE_OBS1 = 'stats_fv3jedi-omb_arrecon' + arrecon_year_str + iop_str + '_exp' + series_id + '0[3456]_bin-' + \
		str(inc_hgt_bins_m) + 'm_nepac_ivt-{below-250,above-250,above-500,above-750}.txt'
infile_name_obs_list = subprocess.check_output('ls ' + PATH_OBS + FILE_OBS1, shell=True, text=True).splitlines()

nfiles = len(infile_name_obs_list)

# Build lists to store required vars
nobs_number_list = []
nobs_percent_list = []
occ_region_name_list = []
occ_dataset_name_list = []
occ_dataset_color_list = []
occ_dataset_linestyle_list = []
occ_stats_type_list = []

# Loop through all the observed occultation files
for ifile in infile_name_obs_list:
	print('')
	print('=======================================================================================')
	print('')

	# Read obs input file
	infile_name_obs = ifile.strip()  # Remove any trailing newlines
	filename_obs = os.path.basename(infile_name_obs)

	# Read in the analysis ASCII file
	names_dtype = ['altitude_bin_mid(km)', 'mean_omb(mrad)', 'mean_omb_percent(%)', 'stddev_omb(mrad)', 'stddev_omb_percent(%)', 
			'rmse_omb(mrad)', 'rmse_omb_percent(%)', 'nobs(#/bin)']
	ncolumns = len(names_dtype)
	formats_dtype = [float] * ncolumns

	# Read the file with explicit field names, using whitespace delimiter
	infile_OBS = np.genfromtxt(
			infile_name_obs,
			dtype=None,
			encoding='utf-8',
			delimiter=None,  # Use whitespace (spaces or tabs) as delimiter
			names=names_dtype,  # Explicitly set the field names
			skip_header=1  # Skip the header line
	)

	print('')
	print(' Observational stats read:')
	print('   ' + infile_name_obs)

	# Gather the data we need from ASCII files
	hgt_ell_bin = infile_OBS['altitude_bin_midkm']
	nobs_number = infile_OBS['nobsbin']

	# Manipulation of the data
	nobs_max = np.max(nobs_number)
	print(nobs_max)
	nobs_percent = nobs_number / nobs_max * 100

	# Get info from file name
	name_file_split = filename_obs.split("_")
	exp_id_raw = name_file_split[3]
	exp_id = exp_id_raw[-2:]
	ivt_subset_raw = name_file_split[-1]
	ivt_subset = ivt_subset_raw[0:-4]
	print(ivt_subset)

	if ivt_subset == 'ivt-above-250':
		region_str = r'IVT $\geq$ 250'  # Use LaTeX for >= symbol
		occ_dataset_linestyle = 'solid'
	elif ivt_subset == 'ivt-above-500':
		region_str = r'IVT $\geq$ 500'  # Use LaTeX for >= symbol
		occ_dataset_linestyle = 'dashed'
	elif ivt_subset == 'ivt-above-750':
		region_str = r'IVT $\geq$ 750'  # Use LaTeX for >= symbol
		occ_dataset_linestyle = 'dotted'
	elif ivt_subset == 'ivt-below-250':
		region_str = r'IVT $<$ 250'
		occ_dataset_linestyle = 'dashdot'

	if exp_id == '03':
		occ_dataset_name = 'ROPP1D'
		occ_dataset_color = 'purple'
	elif exp_id == '04':
		occ_dataset_name = 'ROPP2D'
		occ_dataset_color = 'blue'
	elif exp_id == '05':
		occ_dataset_name = 'NBAM-E'
		occ_dataset_color = 'orange'
	elif exp_id == '06':
		occ_dataset_name = 'NBAM-O'
		occ_dataset_color = 'red'
	else:
		print('')
		print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
		print(' Experiment ID not recognized:')
		print('   ' + exp_id)
		print(' Therefore ABORTING now!!')
		print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
		print('')
		sys.exit()

	# Append to output lists
	nobs_number_list.append(nobs_number)
	nobs_percent_list.append(nobs_percent)
	occ_region_name_list.append(region_str)
	occ_dataset_name_list.append(occ_dataset_name)
	occ_dataset_color_list.append(occ_dataset_color)
	occ_dataset_linestyle_list.append(occ_dataset_linestyle)

print(occ_dataset_name_list)
print(occ_dataset_color_list)
print(occ_dataset_linestyle_list)
print(occ_region_name_list)
print(nobs_number_list[2])
print(nobs_percent_list[2])

# Plotting
# Set up the plot
f, ax = plt.subplots(1, sharey=True, gridspec_kw={'width_ratios': [1]})

linewidth = 1.0

# Plot the data
for ifile in range(nfiles):
	label_str_nobs = f"{occ_region_name_list[ifile]} {occ_dataset_name_list[ifile]}"
	#label_str_nobs = occ_dataset_name_list[ifile] + ' ' + occ_region_name_list[ifile]
	ax.plot(nobs_percent_list[ifile], hgt_ell_bin, linestyle=occ_dataset_linestyle_list[ifile], lw=linewidth,
			label=label_str_nobs, color=occ_dataset_color_list[ifile])

# Define plot axis etc.
top_plot_hgt_km = 8.1
ax.grid(visible=True, which='both', linestyle=':', color='lightgrey')
ax.axis((0, 102, -0.1, top_plot_hgt_km))
ax.yaxis.set_major_locator(plt.MultipleLocator(1))
ax.xaxis.set_major_locator(plt.MultipleLocator(10))
ax.tick_params(axis='both', which='major', labelsize=8)

# Bring subplots close to each other
f.subplots_adjust(hspace=0.0)

ax.set_ylabel('Mean Sea Level Altitude  [km]', fontsize=10)
#ax.set_ylabel('Geometric Height  [km]', fontsize=11)
ax.set_xlabel('Proportion of Total Occultations  [% bin$^{-1}$]', fontsize=10)

ratio = 1.0
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)

# Add a legend
ax.legend(loc='upper left', prop={'size': 8})

# Make everything beautiful automatically
plt.tight_layout(rect=[0, 0.03, 1, 0.90])

# Save the plot
if not os.path.exists(output_path_pdf):
	os.makedirs(output_path_pdf)

output_name_pdf = 'stats_compare_fv3jedi-omb_compare-count-propor-ivt-' + str(inc_hgt_bins_m) + 'm.pdf'
plt.savefig(output_path_pdf + output_name_pdf)

print('')
print(' Plot created:')
print('   ' + output_path_pdf + output_name_pdf)
print('')


