#________________________________________________________________
#
#   Python script to read JEDI hofx output in its native
#    netCDF format and calculate statistics on RO obs. 
#    This example also reads in a IVT derived from the 
#    ERA5 reanalysis in ascii format and subsets the
#    stats based on a given threshold of IVT.
#
#   ~ mjm 2023-02-08
#________________________________________________________________

#import h5py
#import netCDF4
import numpy as np
np.set_printoptions(threshold=np.inf) # make it print out the full array
import matplotlib
matplotlib.use('pdf') # reuired for interactive backend node
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
#import wrf
#import metpy.calc as mcalc
#from metpy.units import units
import os
import sys
import subprocess
#sys.path.insert(0, './geoid_correct')
#from geoid import GeoidHeight
from datetime import datetime, timedelta



 # Read args
nargs = len(sys.argv) - 1
script_name = sys.argv[0] # this will always work.

nargs_req = 6
if ( nargs == nargs_req ):
	arrecon_year_str = sys.argv[1]
	iop_str = sys.argv[2]
	series_id = sys.argv[3]
	exp_id = sys.argv[4]
	region = sys.argv[5]
	io_root = sys.argv[6]
else:
	print("")
	print("")
	print("   Wrong number of arguements, ABORTING!!")
	print("")
	print("   You gave " +str(nargs)+ " of the "+str(nargs_req)+" required args" )
	print("")
	print("     Usage:")
	print("       " + script_name + " <arrecon_year> <iop_str> <series_id> <exp_id> <region> <io_root>")
	print("")
	print("     Example:")
	print("       python " + script_name + " 2022 full 01 02 nepac ../")
	print("")
	#print("     NOTE: the stats arg can be")
	#print("            1) full, 2) inside-ar, or 3) outside-ar")
	print("")
	print("")

	sys.exit()



 # User input
#io_root = ".."
#input_root = '/cw3e/mead/projects/cwp106/projects/hiaper/mjmurphy/modeling/output_fv3jedi/'
year_str = arrecon_year_str


PATH_LIST = io_root + '/stats_ascii/'
FILE_LIST = 'list_occs_arrecon'+year_str+iop_str+'_'+region+'_exp'+series_id+exp_id+'.txt'
infile_name_list = subprocess.check_output('ls '+PATH_LIST+FILE_LIST, shell=True).splitlines()
#infile_name_list = subprocess.check_output('ls '+PATH_LIST+FILE_LIST, shell=True).strip()
#print(infile_name_list)
infile_name_list = [ x.decode('utf-8') for x in infile_name_list ] # needed with python 3
infile_name_list = infile_name_list[0]
#print(infile_name_list)

 # Read in the list ascii file
nlines_header = 1
names_dtype = [ 'occ_id', 'type_id', 'year', 'month', 'day', 'hour', 'minute', 'second', 'lat', 'lon', 'height_km', 'year_round', 'month_round', 'day_round', 'hour_round', 'lat_round', 'lon_round', 'ivt(kg/m/s)', 'iwv(mm)']

ncolumns = len(names_dtype)
#formats_dytpe = [np.float] * ncolumns 	# build list with the same element repeated for the given number of times
formats_dytpe = [ 'U4', np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float ] 

dtypeOBS = {      'names': (names_dtype),
	                'formats': (formats_dytpe) }

#infile_OBS  = np.genfromtxt(infile_name_list, dtype = None, delimiter=",", skip_header=nlines_header,  names=names_dtype)
infile_OBS  = np.loadtxt(infile_name_list, dtype = dtypeOBS, skiprows=nlines_header, delimiter=",", usecols=range(ncolumns))

print('')
print(' List file read:')
print('   ' + infile_name_list)

 # Gather the data we need from ascii files
occ_id_list = infile_OBS['occ_id']
type_id_list = infile_OBS['type_id']
year_list   = infile_OBS['year']
month_list  = infile_OBS['month']
day_list    = infile_OBS['day']
hour_list   = infile_OBS['hour']
minute_list = infile_OBS['minute']
second_list = infile_OBS['second']
ivt_list   = infile_OBS['ivt(kg/m/s)']
iwv_list   = infile_OBS['iwv(mm)']

#print(type_id_list)
#sys.exit()


values, counts = np.unique(type_id_list, return_counts=True)

noccs_total = len(type_id_list)

print('')
#print('----------------------------------------------------------------------------')
print('-------------')
print(values)
print(counts)
print('Total Count = '+ str(noccs_total))


idx_ivt = np.where( (ivt_list < 250) )
ivt_below_250 = ivt_list[idx_ivt]
count_ivt_below_250 = ivt_below_250.shape
count_ivt_below_250 = count_ivt_below_250[0]

print('--------------')
print(' Count IVT below ' + str(250))
print(count_ivt_below_250)

for iivt in range(250,1000,250):
	print('--------------')
	print(' Count IVT above ' + str(iivt))
	idx_ivt = np.where( (ivt_list >= iivt) )
	ivt_above = ivt_list[idx_ivt]
	count_ivt_above = ivt_above.shape
	count_ivt_above = count_ivt_above[0]
	print(count_ivt_above)



sys.exit()

values_cosmic = values[idx_cosmic]
counts_cosmic = counts[idx_cosmic]
noccs_cosmic = np.sum(counts_cosmic)
propor_cosmic = (noccs_cosmic / float(noccs_total)) * 100
print('')
print('-------------')
print('Cosmic2')
print('-------------')
print(values_cosmic)
print(counts_cosmic)
print('Count = '+ str(noccs_cosmic))
print('Proportion of Total = '+ str(propor_cosmic)+' %')

idx_noncos = np.where( (values < 750) | (values > 755 ) )
#print(idx_noncos)
values_noncos = values[idx_noncos]
counts_noncos = counts[idx_noncos]
noccs_noncos = np.sum(counts_noncos)
propor_noncos = (noccs_noncos / float(noccs_total)) * 100
print('')
print('-------------')
print('All other constellations')
print('-------------')
print(values_noncos)
print(counts_noncos)
print('Count = '+ str(noccs_noncos))
print('Proportion of Total = '+ str(propor_noncos)+' %')
print('')



