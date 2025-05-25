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

import h5py
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

nargs_req = 8
if ( nargs == nargs_req ):
	arrecon_year_str = sys.argv[1]
	iop_str = sys.argv[2]
	series_id = sys.argv[3]
	exp_id = sys.argv[4]
	region = sys.argv[5]
	threshold_var = sys.argv[6]
	threshold_value = sys.argv[7]
	#threshold_type = sys.argv[8]
	output_root = sys.argv[8]
else:
	print("")
	print("")
	print("   Wrong number of arguements, ABORTING!!")
	print("")
	print("   You gave " +str(nargs)+ " of the "+str(nargs_req)+" required args" )
	print("")
	print("     Usage:")
	print("       " + script_name + " <arrecon_year> <iop_str> <series_id> <exp_id> <region> <threshold_var> <threshold_value> <output_root>")
	#print("       " + script_name + " <arrecon_year> <iop_str> <series_id> <exp_id> <region> <threshold_var> <threshold_value> <threshold_type> <output_root>")
	print("")
	print("     Example:")
	print("       python " + script_name + " 2022 full 05 04 nepac ivt 250 ../")
	print("")
	print("")
	print("")

	sys.exit()




 # User input
input_root = '/discover/nobackup/mjmurph4'
input_root_path = input_root + '/projects/hiaper/mjmurphy/modeling/output_fv3jedi/'
#aoutput_root = '../'
#iop_num = int(iop_str[-2:]) # last 2 characters are iop  number
#threshold_top = 15 # in km, cuttoff for top of data
threshold_top = 30 # in km, cuttoff for top of data
inc_hgt_bins_m = 200
#inc_hgt_bins_m = 500
#threshold_value = 750.
#threshold_value = 500.
#threshold_value = 250.

#extent = [-179, -113, 15, 60] # AR Recon
#lat_max_thres = 60
#lat_max_thres = 55
#lat_min_thres = 20
#lon_min_thres = -179
#lon_max_thres = -125


 # Set lat/lon extent based on region
if ( region == "tropics"):
	lat_max_thres = 20
	lat_min_thres = -20
	lon_min_thres = -180
	lon_max_thres = 180

elif ( region == "nhml"):
	lat_max_thres = 55
	lat_min_thres = 35
	lon_min_thres = -180
	lon_max_thres = 180

elif ( region == "shml"):
	lat_max_thres = -35
	lat_min_thres = -55
	lon_min_thres = -180
	lon_max_thres = 180

elif ( region == "nepac"):
	lat_max_thres = 55
	lat_min_thres = 20
	lon_min_thres = -179
	lon_max_thres = -125

else:
	print("")
	print(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
	print("   Region not recognized; therefore,")
	print("    ABORTING now!!")
	print("   "+region)
	print(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
	print("")
	sys.exit()


 # User defined function
def round_nearest(x, base=25):
    return base * round(float(x) / base)

def hour_rounder(t):
    # Rounds to nearest hour by adding a timedelta hour if minute >= 30
    return (t.replace(second=0, microsecond=0, minute=0, hour=t.hour)
               +timedelta(hours=t.minute//30))


 # Read in input data 
threshold_value = float(threshold_value) # required because its fed in as an arg

year_str = arrecon_year_str
#PATH_OBS = input_root_path + '/'+year_str+iop_str+'/series'+series_id+'/'+exp_id+'/'+year_str+'-01-11_????/hofx/'
#
#PATH_OBS = input_root_path + '/'+year_str+iop_str+'/series'+series_id+'/'+exp_id+'/'+year_str+'-01-1?_????/hofx/'
PATH_OBS = input_root_path + '/'+year_str+iop_str+'/series'+series_id+'/'+exp_id+'/'+year_str+'-??-??_????/hofx/'
FILE_OBS = 'hofx_gnssro-gfs_'+year_str+'????_????Z.nc4'
infile_name_obs_list = subprocess.check_output('ls '+PATH_OBS+FILE_OBS, shell=True).splitlines()
#print(infile_name_obs_list)
infile_name_obs_list = [ x.decode('utf-8') for x in infile_name_obs_list ] # needed with python 3
#print(infile_name_obs_list)
#sys.exit()

PATH_LIST = output_root + 'stats_ascii/'
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

#print(type(year_list[0]))
noccs_list_file = len(occ_id_list)
#print(noccs_list_file)
date_fmt_list = []
for idate in range(noccs_list_file):
	year_hold_str  = str(int(year_list[idate])).zfill(4)
	month_hold_str = str(int(month_list[idate])).zfill(2)
	day_hold_str   = str(int(day_list[idate])).zfill(2)
	hour_hold_str  = str(int(hour_list[idate])).zfill(2)
	minute_hold_str = str(int(minute_list[idate])).zfill(2)
	second_hold_str = str(int(second_list[idate])).zfill(2)
	date_fmt_hold = year_hold_str +'-'+ month_hold_str +'-'+ day_hold_str +' '+ \
			hour_hold_str +':'+ minute_hold_str +':'+ second_hold_str
	date_fmt_list.append(date_fmt_hold)
date_fmt_list = np.array(date_fmt_list) # conver to numpy array
#print(date_fmt_list)

 # -- Put the data into height bins
threshold_top_m = int(threshold_top * 1000)
hgt_bins_lower_m = range(0,threshold_top_m,inc_hgt_bins_m)
#hgt_bins_mid_m   = range(0+(inc_hgt_bins_m/2),20000+(inc_hgt_bins_m/2),inc_hgt_bins_m)
#hgt_bins_upper_m = range(0+inc_hgt_bins_m,20000+inc_hgt_bins_m,inc_hgt_bins_m)
hgt_bins_upper_m = range(0+inc_hgt_bins_m,threshold_top_m+inc_hgt_bins_m,inc_hgt_bins_m)
#hgt_bins_upper_m = hgt_bins_lower_m + inc_hgt_bins_m

nbins = len(hgt_bins_lower_m)
#print(type(inc_hgt_bins_m))
#print(inc_hgt_bins_m)
#print(int(inc_hgt_bins_m/2))
hgt_mid_km  = [x/1000. for x in range(0+(int(inc_hgt_bins_m/2)),threshold_top_m+(int(inc_hgt_bins_m/2)),inc_hgt_bins_m)]
#hgt_mid_km  = [x/1000. for x in range(0+(int(inc_hgt_bins_m/2)),20000+(int(inc_hgt_bins_m/2)),inc_hgt_bins_m)]
#hgt_mid_km  = [x/1000. for x in range(0+(inc_hgt_bins_m/2),20000+(inc_hgt_bins_m/2),inc_hgt_bins_m)]
hgt_upper_km_out  = [float(x/1000.) for x in hgt_bins_upper_m ]
hgt_lower_km_out  = [float(x/1000.) for x in hgt_bins_lower_m ]
#print(hgt_mid_km)
#print(nbins)
#sys.exit()


 # Build lists to store required vars (this will be a list of lists!)
hgt_km_occ_above_list = []
ba_omb_above_list = []
ba_omb_percent_above_list = []
ba_err_above_list = []
ba_err_eff_above_list = []
ba_err_percent_above_list = []
ba_omb_binned_above_list = []
ba_omb_percent_binned_above_list = []
ba_err_binned_above_list = []
ba_err_eff_binned_above_list = []
ba_err_percent_binned_above_list = []

hgt_km_occ_below_list = []
ba_omb_below_list = []
ba_omb_percent_below_list = []
ba_err_below_list = []
ba_err_eff_below_list = []
ba_err_percent_below_list = []
ba_omb_binned_below_list = []
ba_omb_percent_binned_below_list = []
ba_err_binned_below_list = []
ba_err_eff_binned_below_list = []
ba_err_percent_binned_below_list = []


 # Loop through all the hofx files to get occultations
for ifile in infile_name_obs_list:
	print('')
	print('=======================================================================================')

	infile_name_obs = ifile
	print(infile_name_obs)
	infile_OBS = h5py.File(infile_name_obs,'r')
	#infile_OBS = netCDF4.Dataset(infile_name_obs,'r')
	print('')
	print(' Observation file:')
	print('   ' + infile_name_obs)


     # Gather the data we need
	ba_hofx_raw = infile_OBS['hofx']['bendingAngle'][:]
	ba_obs_raw  = infile_OBS['ObsValue']['bendingAngle'][:]

	date_time_raw = infile_OBS['MetaData']['dateTime'][:] # date in seconds since 1970-01-01T00:00:00Z
	hgt_obs_raw = infile_OBS['MetaData']['height'][:] # Geometric altitude in m
	lat_obs_raw = infile_OBS['MetaData']['latitude'][:] 
	lon_obs_raw = infile_OBS['MetaData']['longitude'][:] # west is negative
	lon_obs_raw = infile_OBS['MetaData']['longitude'][:] # west is negative
	id_sat_raw  = infile_OBS['MetaData']['satelliteIdentifier'][:] # mission id
	id_ref_raw  = infile_OBS['MetaData']['satelliteTransmitterId'][:] # PRN number
	flag_rising_raw = infile_OBS['MetaData']['satelliteAscendingFlag'][:] # rising/setting
	flag_const_raw  = infile_OBS['MetaData']['satelliteConstellationRO'][:] # GPS=401, GLONASS=402, Gal=???
	rec_num_raw = infile_OBS['MetaData']['sequenceNumber'][:] 

      # -- Get metadata
	ba_obs_attrs  = infile_OBS['ObsValue']['bendingAngle'].attrs
	ba_hofx_attrs = infile_OBS['hofx']['bendingAngle'].attrs
	ba_obs_missing_value  = ba_obs_attrs['_FillValue'][0]
	ba_hofx_missing_value = ba_hofx_attrs['_FillValue'][0]

	# Must convert to float128 for consistency
	ba_obs_missing_value  = np.float128(ba_obs_missing_value)
	ba_hofx_missing_value = np.float128(ba_hofx_missing_value)



     # Manipulation of the data
	index_hgt_raw = np.where( (hgt_obs_raw <= threshold_top_m) )
	#index_hgt_raw = np.where( hgt_obs_raw <= 20000)
	index_hgt = index_hgt_raw[0]
	#print(index_hgt)

     # Subset in height
	lat_obs_tmp = lat_obs_raw[index_hgt]
	lon_obs_tmp = lon_obs_raw[index_hgt]
	hgt_obs_tmp = hgt_obs_raw[index_hgt]
	rec_num_tmp = rec_num_raw[index_hgt]
	date_time_tmp = date_time_raw[index_hgt]
	id_sat_tmp = id_sat_raw[index_hgt]
	id_ref_tmp = id_ref_raw[index_hgt]
	flag_rising_tmp = flag_rising_raw[index_hgt]
	flag_const_tmp  = flag_const_raw[index_hgt]

	ba_obs_tmp  = ba_obs_raw[index_hgt]
	ba_hofx_tmp = ba_hofx_raw[index_hgt]
	
     # Get unique occ identifier
	unique_id_occ = np.unique(rec_num_tmp) # record number seems to do the trick
	print(unique_id_occ)
	
	rec_num = rec_num_tmp
	lat_obs = lat_obs_tmp
	lon_obs = lon_obs_tmp
	hgt_obs = hgt_obs_tmp
	ba_obs  = np.float128(ba_obs_tmp)
	#ba_err  = ba_err_tmp
	#ba_err_eff  = ba_err_eff_tmp
	ba_hofx = np.float128(ba_hofx_tmp)

	date_time_obs = date_time_tmp
	id_sat_obs = id_sat_tmp
	id_ref_obs = id_ref_tmp
	flag_rising_obs = flag_rising_tmp
	flag_const_obs  = flag_const_tmp


     # Mask the missing values
	#ba_hofx_missing_value = -3.3687953e+38 # discovered by trial and error
	ba_hofx_mask = np.ma.masked_where( ba_hofx == ba_hofx_missing_value, ba_hofx )
	#sys.exit()

     # Calculate the main variables
	ba_omb = (ba_obs - ba_hofx_mask ) * 1000 # convert to mrad
	ba_omb_percent = ( (ba_obs - ba_hofx_mask) / ba_hofx_mask ) * 100
	hgt_obs_km  = hgt_obs / 1000

     # Save the occultations we need 
	for iid in unique_id_occ:
		lat_tmp = lat_obs[rec_num == iid]
		lon_tmp = lon_obs[rec_num == iid]
		hgt_tmp = hgt_obs[rec_num == iid]

		lat_top = lat_tmp[-1]
		lat_bot = lat_tmp[0]
		lon_top = lon_tmp[-1]
		lon_bot = lon_tmp[0]
		#hgt_top = hgt_tmp[-1]
		hgt_bot = hgt_tmp[0]


	     # Need occs within our region, checking bottom of occ
		if ( ((lat_bot < lat_max_thres) and (lat_bot > lat_min_thres)) \
			and  ((lon_bot < lon_max_thres) and (lon_bot > lon_min_thres)) \
			):
		#if ( ( ((lat_top < lat_max_thres) and (lat_top > lat_min_thres)) \
		#	or ((lat_bot < lat_max_thres) and (lat_bot > lat_min_thres)) ) \
		#	and  ( ((lon_top < lon_max_thres) and (lon_top > lon_min_thres)) \
		#	or ((lon_bot < lon_max_thres) and (lon_bot > lon_min_thres)) ) \
		#	):
			print('------------')
			print('Occ Record Number: ' + str(iid))
			print('bottom at ' + str(hgt_bot) + ' m')
			print(lat_bot)
			print(lon_bot)
			#print(lat_top)
			#print(lon_top)

			hgt_obs_km_sel = hgt_obs_km[rec_num == iid]
			ba_omb_sel = ba_omb[rec_num == iid]
			ba_omb_percent_sel = ba_omb_percent[rec_num == iid]

			date_time_sel = date_time_obs[rec_num == iid]
			id_sat_sel = id_sat_obs[rec_num == iid]
			id_ref_sel = id_ref_obs[rec_num == iid]
			flag_rising_sel = flag_rising_obs[rec_num == iid]
			flag_const_sel = flag_const_obs[rec_num == iid]

			hgt_obs_km_bot = hgt_obs_km_sel[0]
			date_time_bot = date_time_sel[0] # NOTE: time is the same for all hieghts!
			id_sat_bot = id_sat_sel[0] # NOTE: time is the same for all hieghts!
			id_ref_bot = id_ref_sel[0] # NOTE: time is the same for all hieghts!
			flag_rising_bot = flag_rising_sel[0] # NOTE: time is the same for all hieghts!
			flag_const_bot = flag_const_sel[0] # NOTE: time is the same for all hieghts!

		     # Build satellite identifier
			if (flag_rising_bot == 0):
				rising_symbol = "r"
			elif (flag_rising_bot == 1):
				rising_symbol = "s"
			else:
				print('')
				print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
				print(' The flag for rising/setting is not recognized"')
				print(flag_rising_bot)
				print(' therefore, ABORTING now!!')
				print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
				print('')
				sys.exit()


			if (flag_const_bot == 401):
				const_type = "g"
			elif (flag_const_bot == 402):
				const_type = "r"
			else:
				print('')
				print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
				print(' The flag for constellation is not recognized"')
				print(flag_const_bot)
				print(' therefore, ABORTING now!!')
				print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
				print('')
				sys.exit()

			type_dataset = id_sat_bot

			id_ref_bot_str = "%02d" % (id_ref_bot,)
			occ_id = const_type + str(id_ref_bot_str) + rising_symbol
			print('Occultation ID = ' + occ_id)

		     # Time manipulation
			#junk = datetime.utcfromtimestamp(1674832653).strftime('%Y-%m-%d %H:%M:%S')
			date_occ_full = datetime.utcfromtimestamp(date_time_bot).strftime('%Y-%m-%d %H:%M:%S')
			#date_occ_str = datetime.utcfromtimestamp(date_time_bot).strftime('%Y-%m-%d')
			#time_occ_str = datetime.utcfromtimestamp(date_time_bot).strftime('%H%M')
			#date_occ_fmt = time_occ_str + " UTC " + date_occ_str
			occ_year = datetime.utcfromtimestamp(date_time_bot).strftime('%Y')
			occ_month = datetime.utcfromtimestamp(date_time_bot).strftime('%m')
			occ_day = datetime.utcfromtimestamp(date_time_bot).strftime('%d')
			occ_hour = datetime.utcfromtimestamp(date_time_bot).strftime('%H')
			occ_minute = datetime.utcfromtimestamp(date_time_bot).strftime('%M')
			occ_second = datetime.utcfromtimestamp(date_time_bot).strftime('%S')

			idx_occ_list = np.where( (occ_id_list == occ_id) & (date_fmt_list == date_occ_full) )
			print(date_occ_full)
			#print(occ_id)
			#print(type_dataset)

			#print(idx_occ_list)
			occ_ivt = ivt_list[idx_occ_list][0]
			print('IVT = '+ str(occ_ivt))
			occ_iwv = iwv_list[idx_occ_list][0]
			print('IWV = '+ str(occ_iwv))

			if (threshold_var == 'ivt' ):
				occ_subset = occ_ivt
			elif (threshold_var == 'iwv' ):
				occ_subset = occ_iwv
			else:
				print('')
				print(' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
				print('   Chosen threshold variable not recognized,')
				print('    therefore, ABORTING now!!')
				print('   ' + threshold_var)
				print(' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
				print('')
				sys.exit()



		      # Bin the data for the given occultation
			#   (taking mean of values in the given occ over the height range of the given bin)
			ba_omb_binned = np.ma.masked_all(shape=(nbins),dtype='float128') # initialize a fully masked array
			ba_omb_percent_binned = np.ma.masked_all(shape=(nbins),dtype='float128') # initialize a fully masked array
			#ba_omb_binned = np.ma.masked_all(shape=(nbins)) # initialize a fully masked array
			#ba_omb_percent_binned = np.ma.masked_all(shape=(nbins)) # initialize a fully masked array
			for ibin in range(nbins):
				#print('------------')
				#print(ibin)
				hgt_lower_km = float(hgt_bins_lower_m[ibin]) / 1000.
				hgt_upper_km = float(hgt_bins_upper_m[ibin]) / 1000.

				idx_hgt_bin_occ = np.where( (hgt_obs_km_sel > hgt_lower_km) & (hgt_obs_km_sel < hgt_upper_km) )

			     # Test for case with empty height subset
				size_idx = idx_hgt_bin_occ[0].size
				if ( size_idx == 0 ):
					continue

				ba_omb_tmp = ba_omb_sel[idx_hgt_bin_occ]
				ba_omb_percent_tmp = ba_omb_percent_sel[idx_hgt_bin_occ]
		
				ndims_tmp = ba_omb_percent_tmp.shape[0]
				if (  ndims_tmp == 0 ): # just stop here if its empty
					del(ba_omb_percent_tmp)
					continue

				ba_omb_binned[ibin] = np.ma.mean(ba_omb_tmp, axis=0)
				ba_omb_percent_binned[ibin] = np.ma.mean(ba_omb_percent_tmp, axis=0)

				del(ba_omb_tmp)
				del(ba_omb_percent_tmp)


		      # Save the important vars in lists for the given occultation
			if ( occ_subset >= threshold_value):
				print('IVT >= '+ str(threshold_value))

				hgt_km_occ_above_list.append( hgt_obs_km[rec_num == iid] )
				ba_omb_above_list.append( ba_omb[rec_num == iid] )
				ba_omb_percent_above_list.append( ba_omb_percent[rec_num == iid] )

				ba_omb_binned_above_list.append(ba_omb_binned)
				ba_omb_percent_binned_above_list.append(ba_omb_percent_binned)

			else:
				print('IVT < '+ str(threshold_value))

				hgt_km_occ_below_list.append( hgt_obs_km[rec_num == iid] )
				ba_omb_below_list.append( ba_omb[rec_num == iid] )
				ba_omb_percent_below_list.append( ba_omb_percent[rec_num == iid] )

				ba_omb_binned_below_list.append(ba_omb_binned)
				ba_omb_percent_binned_below_list.append(ba_omb_percent_binned)



noccs_above = len(hgt_km_occ_above_list)
noccs_below = len(hgt_km_occ_below_list)
print('')
print('================================')
print('   Total number of occs above threshold = ' + str(noccs_above))
print('   Total number of occs below threshold = ' + str(noccs_below))
#print(hgt_km_occ_list)


 # -- Count values
count_missing_bins_above = np.ma.zeros(shape=(nbins)) # initialize an array full of zeros
for iocc in range(noccs_above):
	for ibin in range(nbins):
		count_missing_bins_above[ibin] = count_missing_bins_above[ibin] + \
						np.ma.count_masked(ba_omb_binned_above_list[iocc][ibin], axis=None)
						#np.ma.count_masked(refrac_diff_binned_list[iocc][ibin], axis=0)
count_obs_bins_above = noccs_above - count_missing_bins_above
#print("Number of occs per bin above threshold:")
#print(count_obs_bins_above)

count_missing_bins_below = np.ma.zeros(shape=(nbins)) # initialize an array full of zeros
for iocc in range(noccs_below):
	print(iocc)
	for ibin in range(nbins):
		count_missing_bins_below[ibin] = count_missing_bins_below[ibin] + \
						np.ma.count_masked(ba_omb_binned_below_list[iocc][ibin], axis=None)
						#np.ma.count_masked(refrac_diff_binned_list[iocc][ibin], axis=0)
count_obs_bins_below = noccs_below - count_missing_bins_below
#print("Number of occs per bin below threshold:")
#print(count_obs_bins_below)
print('================================')
print('')


 # Square of the OmB values for calculation of RMSE later
ba_omb_binned_above_squared_list = []
ba_omb_percent_binned_above_squared_list = []
#ba_omb_binned = np.ma.masked_all(shape=(nbins)) # initialize a fully masked array
for iocc in range(noccs_above):
	ba_omb_binned_above_squared_list.append(np.ma.power(ba_omb_binned_above_list[iocc],2))
	ba_omb_percent_binned_above_squared_list.append(np.ma.power(ba_omb_percent_binned_above_list[iocc],2))

ba_omb_binned_below_squared_list = []
ba_omb_percent_binned_below_squared_list = []
for iocc in range(noccs_below):
	ba_omb_binned_below_squared_list.append(np.ma.power(ba_omb_binned_below_list[iocc],2))
	ba_omb_percent_binned_below_squared_list.append(np.ma.power(ba_omb_percent_binned_below_list[iocc],2))


 # -- Calculate the mean stats
ba_omb_above_mean = np.ma.mean(ba_omb_binned_above_list, axis=0)
ba_omb_above_std = np.ma.std(ba_omb_binned_above_list, axis=0)
ba_omb_above_rmse = np.ma.sqrt(np.ma.mean(ba_omb_binned_above_squared_list, axis=0))
#ba_omb_above_rmse = np.ma.sqrt(np.ma.mean(ba_omb_binned_above_list, axis=0)**2)
ba_omb_percent_above_mean = np.ma.mean(ba_omb_percent_binned_above_list, axis=0)
ba_omb_percent_above_std = np.ma.std(ba_omb_percent_binned_above_list, axis=0)
ba_omb_percent_above_rmse = np.ma.sqrt(np.ma.mean(ba_omb_percent_binned_above_squared_list, axis=0))
#ba_omb_percent_above_rmse = np.ma.sqrt(np.ma.mean(ba_omb_percent_binned_above_list, axis=0)**2)

ba_omb_below_mean = np.ma.mean(ba_omb_binned_below_list, axis=0)
ba_omb_below_std = np.ma.std(ba_omb_binned_below_list, axis=0)
ba_omb_below_rmse = np.ma.sqrt(np.ma.mean(ba_omb_binned_below_squared_list, axis=0))
#ba_omb_below_rmse = np.ma.sqrt(np.ma.mean(ba_omb_binned_below_list, axis=0)**2)
ba_omb_percent_below_mean = np.ma.mean(ba_omb_percent_binned_below_list, axis=0)
ba_omb_percent_below_std = np.ma.std(ba_omb_percent_binned_below_list, axis=0)
ba_omb_percent_below_rmse = np.ma.sqrt(np.ma.mean(ba_omb_percent_binned_below_squared_list, axis=0))
#ba_omb_percent_below_rmse = np.ma.sqrt(np.ma.mean(ba_omb_percent_binned_below_list, axis=0)**2)

#print("mean")
#print(hgt_mid_km)
#print(ba_omb_percent_mean)
#print(ba_omb_percent_std)
#print(ba_omb_percent_rmse)
#print(len(ba_omb_mean))
#print(len(ba_omb_rmse))
#print(len(ba_omb_percent_rmse))
#sys.exit()

# Plotting
 #-- Set up the plot
#f, axarr = plt.subplots()
#f, axarr = plt.subplots(1, 2, sharey=True, gridspec_kw = {'width_ratios':[1, 0.25]})
f, axarr = plt.subplots(1, 2, sharey=True, gridspec_kw = {'width_ratios':[1, 0.3]})
#f, axarr = plt.subplots(1, 3, sharey=True, gridspec_kw = {'width_ratios':[1, 0.3, 0.3]})
#title_plot = 'AR Recon All Years All IOPs Dropsondes\n' + \
#		stats_type.capitalize() + ' Statistics versus Reanalysis'
#title_plot = 'AR Recon '+arrecon_year_str+' '+iop_str.upper()+' Series '+series_id+'.'+exp_id + '\n' + \
title_plot = 'Experiment '+series_id+"."+exp_id.upper() + ' AR Recon ' + year_str + ' '+ region.capitalize() +' '+ \
		threshold_var.upper() + ' >=' + str(threshold_value)
#title_plot = 'Experiment ' +exp_id.upper() + ' at ' + hour_str + ' UTC ' + date_str
f.suptitle(title_plot)

 #-- Plot the main data
for iocc in range(noccs_above):
	if (iocc == 0):
		axarr[0].plot(ba_omb_percent_above_list[iocc], hgt_km_occ_above_list[iocc], linestyle='-', lw=0.75, label='Occ.', color='darkgrey')
	else:
		axarr[0].plot(ba_omb_percent_above_list[iocc], hgt_km_occ_above_list[iocc], linestyle='-', lw=0.75, label=None, color='darkgrey')

	#print('here')
	#print(noccs)


 #-- Plot reference lines (constant)
#axarr[1].axvline(x=0.0, linestyle='--', color='black')
axarr[0].axvline(x=0.0, linestyle='--', color='black')

 # -- Plot additional stats
color_occ = 'red'
axarr[0].plot(ba_omb_percent_above_mean, hgt_mid_km, linestyle='-', label=r'$\mu$', color=color_occ)
axarr[0].plot(ba_omb_percent_above_std, hgt_mid_km, linestyle=':', label=r'$\sigma$', color=color_occ)

 # -- Plot data on extra panels
axarr[1].plot(count_obs_bins_above, hgt_mid_km, linestyle='-', label='count', color=color_occ)

 # -- Add text stamp with number of occs
stamp_str = 'n = ' + str(noccs_above)
axarr[0].text(16, 16.0, stamp_str, fontsize=9)
#axarr[0].text(60, 16.0, stamp_str, fontsize=9)
#axarr[0].text(60, 12.0, stamp_str, fontsize=9)
#axarr[0].text(210, 11.5, stamp_str, fontsize=10)

 #-- Add a legend
axarr[0].legend(loc='upper right', prop={'size': 8})
axarr[1].legend(loc='upper left', prop={'size': 8})


 #-- Define plot axis etc.
#top_plot_hgt_km = 15.1
top_plot_hgt_km = 20.1
#top_plot_hgt_km = threshold_top

axarr[0].grid(b=True, which='both', linestyle=':', color='lightgrey') # grid lines
axarr[0].axis((-28,28,-0.1,top_plot_hgt_km)) # x1,x2,y1,y2
#axarr[0].axis((-100,100,-0.1,top_plot_hgt_km)) # x1,x2,y1,y2
axarr[0].yaxis.set_major_locator(plt.MultipleLocator(1)) # ticks are a multipule of X
axarr[0].xaxis.set_major_locator(plt.MultipleLocator(5)) # ticks are a multipule of X
#axarr[0].xaxis.set_major_locator(plt.MultipleLocator(20)) # ticks are a multipule of X
axarr[0].tick_params(axis='both', which='major', labelsize=8)

#noccs_min = round( min(count_obs_bins)/10)*10 # round to nearest tenth
import math
noccs_min = math.floor( min(count_obs_bins_above)/10)*10 # round to nearest tenth
print(noccs_min)
axarr[1].grid(b=True, which='both', linestyle=':', color='lightgrey') # grid lines
axarr[1].axis((noccs_min,noccs_above+2,-0.1,top_plot_hgt_km)) # x1,x2,y1,y2
axarr[1].yaxis.set_major_locator(plt.MultipleLocator(1)) # ticks are a multipule of X
#axarr[1].xaxis.set_major_locator(plt.MultipleLocator(10)) # ticks are a multipule of X
axarr[1].tick_params(axis='both', which='major', labelsize=8)

axarr[0].set_ylabel('Altitude  [km]', fontsize=10)
axarr[0].set_xlabel(r'$\alpha$' +' Obs Minus '+ r'$\alpha$' +' Bkgd / '+ r'$\alpha$' +' Bkgd  [%]', fontsize=10)
#axarr.set_xlabel('Obs '+ r'$\alpha$' +' Minus Bkgd '+ r'$\alpha$' +' / Bkgd '+ r'$\alpha$' +' [%]', fontsize=10)
#axarr.set_xlabel('Obs Minus Background '+ r'$\alpha$' +' / Background  [%]', fontsize=10)
axarr[1].set_xlabel('Obs  [#/bin]', fontsize=10)


 #-- Set aspect ratio (larger numbers are more narrow)
ratio = 1.25
#ratio = 2.25


 # Write out to ascii
outfile_header = "altitude_bin_mid(km) mean_omb(mrad) mean_omb_percent(%) stddev_omb(mrad) stddev_omb_percent(%) rmse_omb(mrad) rmse_omb_percent(%) nobs(#/bin)"
output_array = zip(hgt_mid_km, ba_omb_above_mean, ba_omb_percent_above_mean, \
                        ba_omb_above_std, ba_omb_percent_above_std, \
                        ba_omb_above_rmse, ba_omb_percent_above_rmse, \
                        count_obs_bins_above)

output_path_ascii = output_root + 'stats_ascii/'
if not os.path.exists(output_path_ascii):
	os.makedirs(output_path_ascii)

outfile_name_above = output_path_ascii + \
		'stats_fv3jedi-omb_arrecon'+year_str+iop_str+ \
			'_exp'+series_id+exp_id+'_bin-'+str(inc_hgt_bins_m)+'m_'+region+'_'+ \
			threshold_var+'-above-'+str(int(threshold_value))+'.txt'

np.savetxt(outfile_name_above, list(output_array), delimiter=" ", header=outfile_header, fmt='%1.4f')

print('')
print(' Stats file (ascii) created:')
print('   ' + outfile_name_above)
print('')

outfile_header = "altitude_bin_mid(km) mean_omb(mrad) mean_omb_percent(%) stddev_omb(mrad) stddev_omb_percent(%) rmse_omb(mrad) rmse_omb_percent(%) nobs(#/bin)"
output_array = zip(hgt_mid_km, ba_omb_below_mean, ba_omb_percent_below_mean, \
                        ba_omb_below_std, ba_omb_percent_below_std, \
                        ba_omb_below_rmse, ba_omb_percent_below_rmse, \
                        count_obs_bins_below)

output_path_ascii = output_root + 'stats_ascii/'
if not os.path.exists(output_path_ascii):
	os.makedirs(output_path_ascii)

outfile_name_below = output_path_ascii + \
		'stats_fv3jedi-omb_arrecon'+year_str+iop_str+ \
			'_exp'+series_id+exp_id+'_bin-'+str(inc_hgt_bins_m)+'m_'+region+'_'+ \
			threshold_var+'-below-'+str(int(threshold_value))+'.txt'

np.savetxt(outfile_name_below, list(output_array), delimiter=" ", header=outfile_header, fmt='%1.4f')

print('')
print(' Stats file (ascii) created:')
print('   ' + outfile_name_below)
print('')



 # Make the image
output_path_image = output_root + 'plots/'
if not os.path.exists(output_path_image):
	os.makedirs(output_path_image)


output_name_image = 'stats_fv3jedi-omb_arrecon'+year_str+iop_str+ \
			'_exp'+series_id+exp_id+'_bin-'+str(inc_hgt_bins_m)+'m_'+region+'_'+ \
			threshold_var+'-above-'+str(int(threshold_value))+'.pdf'

#output_name_image = 'stats_omb_ro.pdf'
outfile_image = output_path_image + output_name_image
plt.savefig(outfile_image)

print('')
print(' Plot created:')
print('   ' + outfile_image)
print('')


#list_file.close()


sys.exit()


