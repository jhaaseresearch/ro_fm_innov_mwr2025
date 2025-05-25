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
from matplotlib.ticker import StrMethodFormatter

import subprocess
import sys
import os
from datetime import datetime, timedelta
sys.path.insert(0, './calc_vars')
from label_axes_iterate import label_axes

# Read args
nargs = len(sys.argv) - 1
script_name = sys.argv[0]  # this will always work.

nargs_req = 8
if nargs == nargs_req:
	arrecon_year_str = sys.argv[1]
	iop_str = sys.argv[2]
	series_id = sys.argv[3]
	threshold_var = sys.argv[4]
	threshold_value = sys.argv[5]
	threshold_type = sys.argv[6]
	region = sys.argv[7]
	io_root = sys.argv[8]
else:
	print()
	print()
	print("   Wrong number of arguments, ABORTING!!")
	print()
	print(f"   You gave {nargs} of the {nargs_req} required args")
	print()
	print("     Usage:")
	print(f"       {script_name} <arrecon_year> <iop_str> <series_id> <threshold_var> <threshold_value> <threshold_type> <region> <output_root>")
	print()
	print("     Example:")
	print(f"       python {script_name} 2022 full 05 ivt 250 above nepac ../")
	print()
	print()
	sys.exit()

threshold_value = float(threshold_value)

# User input
inc_hgt_bins_m = 200
missingValue = -999.0  # this same missing value seems to be used for all variables in soundings

nlines_header = 1

input_dir_root = io_root
output_dir_root = io_root
output_path_pdf = output_dir_root + 'plots/'

# Make list of input files
PATH_OBS = input_dir_root + 'stats_ascii/'
FILE_OBS = 'stats_fv3jedi-omb_arrecon' + arrecon_year_str + iop_str + '_exp' + series_id + '0[3456]_bin-' + str(inc_hgt_bins_m) + 'm_' + \
		region + '_' + threshold_var + '-' + threshold_type + '-' + str(int(threshold_value)) + '.txt'

# Use subprocess to list files, decoding output directly
infile_name_obs_list = subprocess.check_output('ls ' + PATH_OBS + FILE_OBS, shell=True, text=True).splitlines()

nfiles = len(infile_name_obs_list)

# Build lists to store required vars
omb_mean_list = []
omb_stddev_list = []
omb_mean_percent_list = []
omb_stddev_percent_list = []
omb_rmse_percent_list = []
nobs_percent_list = []
nobs_number_list = []
occ_dataset_name_list = []
occ_dataset_color_list = []
occ_stats_type_list = []

# Loop through all the observed occultation files
for ifile in infile_name_obs_list:
	print()
	print('=======================================================================================')
	print()

	# Read obs input file
	infile_name_obs = ifile.strip()  # Remove any trailing newlines
	filename_obs = os.path.basename(infile_name_obs)

	# Read in the analysis ASCII file
	names_dtype = ['altitude_bin_mid(km)', 'mean_omb(mrad)', 'mean_omb_percent(%)', 'stddev_omb(mrad)', 'stddev_omb_percent(%)', 'rmse_omb(mrad)', 'rmse_omb_percent(%)', 'nobs(#/bin)']
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

	print()
	print(' Observational stats read:')
	print('   ' + infile_name_obs)

	# Gather the data we need from ASCII files
	# Use the sanitized field names as they are loaded by np.genfromtxt
	omb_mean = infile_OBS['mean_ombmrad']
	omb_stddev = infile_OBS['stddev_ombmrad']
	omb_rmse = infile_OBS['rmse_ombmrad']
	omb_mean_percent = infile_OBS['mean_omb_percent']
	omb_stddev_percent = infile_OBS['stddev_omb_percent']
	omb_rmse_percent = infile_OBS['rmse_omb_percent']
	hgt_ell_bin = infile_OBS['altitude_bin_midkm']
	nobs_number = infile_OBS['nobsbin']

	# Manipulation of the data
	omb_mean = np.ma.asarray(omb_mean)
	omb_mean_percent = np.ma.asarray(omb_mean_percent)
	omb_stddev = np.ma.asarray(omb_stddev)
	omb_stddev_percent = np.ma.asarray(omb_stddev_percent)
	omb_rmse_percent = np.ma.asarray(omb_rmse_percent)

	omb_mean = np.ma.masked_where(omb_mean == missingValue, omb_mean)
	omb_mean_percent = np.ma.masked_where(omb_mean_percent == missingValue, omb_mean_percent)
	omb_stddev = np.ma.masked_where(omb_stddev == missingValue, omb_stddev)
	omb_stddev_percent = np.ma.masked_where(omb_stddev_percent == missingValue, omb_stddev_percent)
	omb_rmse_percent = np.ma.masked_where(omb_rmse_percent == missingValue, omb_rmse_percent)

	nobs_max = np.max(nobs_number)
	nobs_percent = nobs_number / nobs_max * 100

	# Get info from file name
	name_file_split = filename_obs.split("_")
	exp_id_raw = name_file_split[3]
	exp_id = exp_id_raw[-2:]

	if exp_id == '02':
		occ_dataset_name = 'NBAM'
		occ_dataset_color = 'red'
	elif exp_id == '03':
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
		print()
		print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
		print(' Experiment ID not recognized:')
		print('   ' + exp_id)
		print(' Therefore ABORTING now!!')
		print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
		print()
		sys.exit()

	print(occ_dataset_name)

	# Append to output lists
	omb_mean_list.append(omb_mean)
	omb_stddev_list.append(omb_stddev)
	omb_mean_percent_list.append(omb_mean_percent)
	omb_stddev_percent_list.append(omb_stddev_percent)
	omb_rmse_percent_list.append(omb_rmse_percent)
	occ_dataset_name_list.append(occ_dataset_name)
	occ_dataset_color_list.append(occ_dataset_color)
	nobs_percent_list.append(nobs_percent)
	nobs_number_list.append(nobs_number)

# Plotting
f, axarr = plt.subplots(1, 3, sharey=True, gridspec_kw={'width_ratios': [1, 1, 1]})

linestyles_mean = ['-', '-', '-', '-', '-', '-', '-', '-']
linestyles_stddev = linestyles_mean
linestyles_nobs = linestyles_mean
linestyles_rmse = linestyles_mean
linewidth = 1.0

# Plot the data
for ifile in range(nfiles):
	label_str_nobs = occ_dataset_name_list[ifile]
	label_str_mean = r'$\mu$'
	label_str_stddev = r'$\sigma$'

	if ifile == 0:
		axarr[0].plot(nobs_number_list[ifile], hgt_ell_bin, linestyle=linestyles_nobs[ifile], lw=linewidth, \
				label=label_str_nobs, color=occ_dataset_color_list[ifile])
		axarr[1].plot(omb_mean_percent_list[ifile], hgt_ell_bin, linestyle=linestyles_mean[ifile], lw=linewidth, \
				label=label_str_mean, color=occ_dataset_color_list[ifile])
		axarr[2].plot(omb_stddev_percent_list[ifile], hgt_ell_bin, linestyle=linestyles_stddev[ifile], lw=linewidth, \
				label=label_str_stddev, color=occ_dataset_color_list[ifile])
	else:
		axarr[0].plot(nobs_number_list[ifile], hgt_ell_bin, linestyle=linestyles_nobs[ifile], lw=linewidth, \
				label=label_str_nobs, color=occ_dataset_color_list[ifile])
		axarr[1].plot(omb_mean_percent_list[ifile], hgt_ell_bin, linestyle=linestyles_mean[ifile], lw=linewidth, \
				label=label_str_mean, color=occ_dataset_color_list[ifile])
		axarr[2].plot(omb_stddev_percent_list[ifile], hgt_ell_bin, linestyle=linestyles_stddev[ifile], lw=linewidth, \
				label=label_str_stddev, color=occ_dataset_color_list[ifile])

# Plot reference lines (constant)
axarr[1].axvline(x=0.0, linestyle='--', color='black')

# Define plot axis etc.
top_plot_hgt_km = 15.1

axarr[0].grid(True, which='both', linestyle=':', color='lightgrey')
axarr[0].yaxis.set_major_locator(plt.MultipleLocator(1))
axarr[0].xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
axarr[0].tick_params(axis='both', which='major', labelsize=7)
#axarr[0].tick_params(axis='both', which='major', labelsize=8)

axarr[1].grid(True, which='both', linestyle=':', color='lightgrey')
axarr[1].axis((-6.7, 6.7, -0.1, top_plot_hgt_km))
axarr[1].yaxis.set_major_locator(plt.MultipleLocator(1))
axarr[1].xaxis.set_major_locator(plt.MultipleLocator(2))
axarr[1].xaxis.set_minor_locator(plt.MultipleLocator(1))
axarr[1].tick_params(axis='both', which='major', labelsize=7)
#axarr[1].tick_params(axis='both', which='major', labelsize=8)

axarr[2].grid(True, which='both', linestyle=':', color='lightgrey')
axarr[2].axis((0, 21.1, -0.1, top_plot_hgt_km))
axarr[2].xaxis.set_major_locator(plt.MultipleLocator(4))
axarr[2].xaxis.set_minor_locator(plt.MultipleLocator(2))
axarr[2].tick_params(axis='both', which='major', labelsize=7)
#axarr[2].tick_params(axis='both', which='major', labelsize=8)

axarr[0].set_ylabel('Mean Sea Level Altitude  [km]', fontsize=8)
root_xlabel_str = r'$\dfrac{O_{\alpha}' + u"\u2212" + r'B_{\alpha}}{B_{\alpha}}$' + '  [%]'
#root_xlabel_str = '(O' + r'$_{\alpha}$' + u"\u2212" + 'B' + r'$_{\alpha}$' + ')' + 'B' + r'$_{\alpha}$' + '$^{-1}$' + '  [%]'
axarr[0].set_xlabel('Count  [tangent pts. bin$^{-1}$]', fontsize=8.0)
axarr[1].set_xlabel('Mean ' + root_xlabel_str, fontsize=8.0)
axarr[2].set_xlabel('Std Dev ' + root_xlabel_str, fontsize=8.0)
#axarr[0].set_xlabel('Count  [tangent pts. bin$^{-1}$]', fontsize=8.5)
#axarr[1].set_xlabel('Mean ' + root_xlabel_str, fontsize=8.5)
#axarr[2].set_xlabel('Std Dev ' + root_xlabel_str, fontsize=8.5)


# Bring subplots close to each other
f.subplots_adjust(wspace=0.1, hspace=0)

# Add a legend
axarr[0].legend(loc='center left', prop={'size': 7})
#axarr[0].legend(loc='center left', prop={'size': 8})

# Add figure panel label strings
label_axes(f, loc=(.26, .94), ha='right', fontsize=13)
#label_axes(f, loc=(.26, .92), ha='right', fontsize=15)

# Save the plot
if not os.path.exists(output_path_pdf):
	os.makedirs(output_path_pdf)

output_name_pdf = 'stats_compare_fv3jedi-omb_count-mean-sigma_bin-' + str(inc_hgt_bins_m) + 'm_' + \
			threshold_var + '-' + threshold_type + '-' + str(int(threshold_value)) + '_' + region + '.pdf'
plt.savefig(output_path_pdf + output_name_pdf)

print('')
print(' Output image created:')
print('    ' + output_path_pdf + output_name_pdf)
print('')

plt.close()

