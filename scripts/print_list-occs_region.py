#________________________________________________________________
#
#   Python script to read JEDI hofx output in its native
#    netCDF format and calculate statistics on RO obs. 
#    This example is set up for plotting the observation
#    minus the background (OmB).
#    
#   ~ mjm 2022-06-29
#________________________________________________________________

import netCDF4
#import h5py
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



 # User input
input_root = '/discover/nobackup/mjmurph4'
input_root_path = input_root + '/projects/hiaper/mjmurphy/modeling/output_fv3jedi/'
io_root = '../'
threshold_top = 30 # in km, cuttoff for top of data
inc_hgt_bins_m = 200
#inc_hgt_bins_m = 500
analysis_name = 'era5'
#analysis_type = 'plevel'


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
	print("       " + script_name + " <arrecon_year> <iop_str> <series_id> <exp_id> <region> <oi_root>")
	print("")
	print("     Example:")
	print("       python " + script_name + " 2022 full 01 02 nepac ../")
	print("")
	#print("     NOTE: the stats arg can be")
	#print("            1) full, 2) inside-ar, or 3) outside-ar")
	print("")
	print("")

	sys.exit()




 # Set lat/lon extent based on region
if ( region == "global"):
	lat_max_thres = 91 # adding 1 to be sure to get edges
	lat_min_thres = -91
	lon_min_thres = -181
	lon_max_thres = 181

elif ( region == "tropics"):
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


# -- strip info from filename
year_str = arrecon_year_str


 # Read in the obs input file
PATH_OBS = input_root_path + '/'+year_str+iop_str+'/series'+series_id+'/'+exp_id+'/'+year_str+'-??-??_????/hofx/'
FILE_OBS = 'hofx_gnssro-gfs_'+year_str+'????_????Z.nc4'
infile_name_obs_list = subprocess.check_output('ls '+PATH_OBS+FILE_OBS, shell=True).splitlines()
#print(infile_name_obs_list)
infile_name_obs_list = [ x.decode('utf-8') for x in infile_name_obs_list ] # needed with python 3
#print(infile_name_obs_list)
#sys.exit()

 # -- Put the data into height bins
threshold_top_m = int(threshold_top * 1000)
hgt_bins_lower_m = range(0,threshold_top_m,inc_hgt_bins_m)
#hgt_bins_mid_m   = range(0+(inc_hgt_bins_m/2),20000+(inc_hgt_bins_m/2),inc_hgt_bins_m)
hgt_bins_upper_m = range(0+inc_hgt_bins_m,20000+inc_hgt_bins_m,inc_hgt_bins_m)
#hgt_bins_upper_m = hgt_bins_lower_m + inc_hgt_bins_m

nbins = len(hgt_bins_lower_m)
#print(type(inc_hgt_bins_m))
#print(inc_hgt_bins_m)
#print(int(inc_hgt_bins_m/2))
hgt_mid_km  = [x/1000. for x in range(0+(int(inc_hgt_bins_m/2)),20000+(int(inc_hgt_bins_m/2)),inc_hgt_bins_m)]
#hgt_mid_km  = [x/1000. for x in range(0+(inc_hgt_bins_m/2),20000+(inc_hgt_bins_m/2),inc_hgt_bins_m)]
hgt_upper_km_out  = [float(x/1000.) for x in hgt_bins_upper_m ]
hgt_lower_km_out  = [float(x/1000.) for x in hgt_bins_lower_m ]
#print(hgt_mid_km)
#print(nbins)
#sys.exit()


 # Build lists to store required vars (this will be a list of lists!)
hgt_km_occ_list = []
ba_omb_list = []
ba_omb_percent_list = []
ba_err_list = []
ba_err_eff_list = []
ba_err_percent_list = []
ba_omb_binned_list = []
ba_omb_percent_binned_list = []
ba_err_binned_list = []
ba_err_eff_binned_list = []
ba_err_percent_binned_list = []



analyses_dir = io_root + 'analyses/'+analysis_name+'/'
output_dir = io_root + 'stats_ascii/'

list_file_name = 'list_occs_arrecon'+arrecon_year_str+iop_str+'_'+region+'_exp'+series_id+exp_id+'.txt'

if not os.path.exists(output_dir):
	os.makedirs(output_dir)

 # Create output file
list_file = f = open(output_dir+list_file_name, 'w')
list_file.write('occ_id, type_id, year, month, day, hour, minute, second, lat, lon, height_km, year_round, month_round, day_round, hour_round, lat_round, lon_round, ivt(kg/m/s), iwv(mm)\n')
#list_file.write('occ_id, type, year, month, day, hour, minute, second, iop, lat, lon, height_km, year_round, month_round, day_round, hour_round, lat_round, lon_round, ivt(kg/m/s), iwv(mm)\n')


 # Loop through all the files to get occultations
for ifile in infile_name_obs_list:
	print('')
	print('=======================================================================================')

	infile_name_obs = ifile
	#infile_OBS = h5py.File(infile_name_obs,'r')
	infile_OBS = netCDF4.Dataset(infile_name_obs,'r')
	print('')
	print(' Observation file:')
	print('   ' + infile_name_obs)

     # Gather the data we need
	infile_group_hofx = infile_OBS['hofx']
	infile_group_obs  = infile_OBS['ObsValue']
	infile_group_err_eff  = infile_OBS['EffectiveError']
	infile_group_err  = infile_OBS['ObsError']
	infile_group_meta = infile_OBS['MetaData']

	hgt_obs_raw = infile_group_meta.variables['height'][:] # Geometric altitude in m
	lat_obs_raw = infile_group_meta.variables['latitude'][:] 
	lon_obs_raw = infile_group_meta.variables['longitude'][:] # west is negative
	rec_num_raw = infile_group_meta.variables['sequenceNumber'][:] # 
	date_time_raw = infile_group_meta.variables['dateTime'][:] # date
	id_sat_raw  = infile_group_meta.variables['satelliteIdentifier'][:] # mission id
	id_ref_raw  = infile_group_meta.variables['satelliteTransmitterId'][:] # PRN number
	flag_rising_raw = infile_group_meta.variables['satelliteAscendingFlag'][:] # rising/setting
	flag_const_raw  = infile_group_meta.variables['satelliteConstellationRO'][:] # GPS=401, GLONASS=402, Gal=???

     # Manipulation of the data
	index_hgt_raw = np.where( (hgt_obs_raw <= threshold_top_m) )
	index_hgt = index_hgt_raw[0]

     # Subset in height
	lat_obs_tmp = lat_obs_raw[index_hgt]
	lon_obs_tmp = lon_obs_raw[index_hgt]
	hgt_obs_tmp = hgt_obs_raw[index_hgt]
	date_time_tmp = date_time_raw[index_hgt]
	rec_num_tmp = rec_num_raw[index_hgt]
	id_sat_tmp = id_sat_raw[index_hgt]
	id_ref_tmp = id_ref_raw[index_hgt]
	flag_rising_tmp = flag_rising_raw[index_hgt]
	flag_const_tmp  = flag_const_raw[index_hgt]

     # Get unique occ identifier
	unique_id_occ = np.unique(rec_num_tmp) # record number seems to do the trick
	#print(unique_id_occ)
	
	rec_num = rec_num_tmp
	lat_obs = lat_obs_tmp
	lon_obs = lon_obs_tmp
	hgt_obs = hgt_obs_tmp
	date_time_obs = date_time_tmp
	id_sat_obs = id_sat_tmp
	id_ref_obs = id_ref_tmp
	flag_rising_obs = flag_rising_tmp
	flag_const_obs  = flag_const_tmp

     # Calculate the main variables
	hgt_obs_km  = hgt_obs / 1000

     # Save the occultations we need 
	for iid in unique_id_occ:
		lat_tmp = lat_obs[rec_num == iid]
		lon_tmp = lon_obs[rec_num == iid]
		hgt_tmp = hgt_obs[rec_num == iid]

		id_sat_sel = id_sat_obs[rec_num == iid]
		id_ref_sel = id_ref_obs[rec_num == iid]
		flag_rising_sel = flag_rising_obs[rec_num == iid]
		flag_const_sel = flag_const_obs[rec_num == iid]
		date_time_sel = date_time_obs[rec_num == iid]
		hgt_obs_km_sel = hgt_obs_km[rec_num == iid]

		lat_top = lat_tmp[-1]
		lat_bot = lat_tmp[0]
		lon_top = lon_tmp[-1]
		lon_bot = lon_tmp[0]
		hgt_bot = hgt_tmp[0]

		hgt_obs_km_bot = hgt_obs_km_sel[0]
		date_time_bot = date_time_sel[0] # NOTE: time is the same for all hieghts!
		id_sat_bot = id_sat_sel[0] # NOTE: time is the same for all hieghts!
		id_ref_bot = id_ref_sel[0] # NOTE: time is the same for all hieghts!
		flag_rising_bot = flag_rising_sel[0] # NOTE: time is the same for all hieghts!
		flag_const_bot = flag_const_sel[0] # NOTE: time is the same for all hieghts!

	     # Need occs within our region, checking bottom of occ
		if ( ((lat_bot < lat_max_thres) and (lat_bot > lat_min_thres)) \
			and  ((lon_bot < lon_max_thres) and (lon_bot > lon_min_thres)) \
			):
		#if ( ( ((lat_top < lat_max_thres) and (lat_top > lat_min_thres)) \
		#	or ((lat_bot < lat_max_thres) and (lat_bot > lat_min_thres)) ) \
		#	and  ( ((lon_top < lon_max_thres) and (lon_top > lon_min_thres)) \
		#	or ((lon_bot < lon_max_thres) and (lon_bot > lon_min_thres)) ) \
		#	):
			print('')
			print('-------------------------------------')
			print('Occ Record Number: ' + str(iid))
			print('bottom at ' + str(hgt_bot) + ' m')
			print('')

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

			print(str(lat_bot)+'N & '+ str(lon_bot)+'E')
			lat_round_raw = round_nearest(lat_bot*100,base=25)
			lat_round = lat_round_raw/100
			lon_round_raw = round_nearest(lon_bot*100,base=25)
			lon_round = lon_round_raw/100
			print("Rounded to:")
			print(str(lat_round)+'N & '+ str(lon_round)+'E')

		     # Time manipulation
			#junk = datetime.utcfromtimestamp(1674832653).strftime('%Y-%m-%d %H:%M:%S')
			date_occ_full = datetime.utcfromtimestamp(date_time_bot).strftime('%Y-%m-%d %H:%M:%S')
			date_occ_str = datetime.utcfromtimestamp(date_time_bot).strftime('%Y-%m-%d')
			time_occ_str = datetime.utcfromtimestamp(date_time_bot).strftime('%H%M')
			date_occ_fmt = time_occ_str + " UTC " + date_occ_str
			occ_year = datetime.utcfromtimestamp(date_time_bot).strftime('%Y')
			occ_month = datetime.utcfromtimestamp(date_time_bot).strftime('%m')
			occ_day = datetime.utcfromtimestamp(date_time_bot).strftime('%d')
			occ_hour = datetime.utcfromtimestamp(date_time_bot).strftime('%H')
			occ_minute = datetime.utcfromtimestamp(date_time_bot).strftime('%M')
			occ_second = datetime.utcfromtimestamp(date_time_bot).strftime('%S')

			print('')
			print(date_occ_fmt)
			#print(date_occ_full)

			date_occ_obj = datetime.strptime(date_occ_full , '%Y-%m-%d %H:%M:%S')
			date_occ_obj_round = hour_rounder(date_occ_obj)
			#print(date_occ_obj)
			#print(date_occ_obj_round)
			date_occ_analysis = date_occ_obj_round.strftime('%H%M UTC %Y-%m-%d')
			date_occ_analysis_file = date_occ_obj_round.strftime('%Y-%m-%d-%H%M')
			year_round = date_occ_obj_round.strftime('%Y')
			month_round = date_occ_obj_round.strftime('%m')
			day_round = date_occ_obj_round.strftime('%d')
			hour_round = date_occ_obj_round.strftime('%H')
			print("Rounded to:")
			print(date_occ_analysis)

		     # Find reanalysis file
			FILE_ANL = 'era5_sfc_0.25deg-grid_subset-arrecon_'+date_occ_analysis_file+'.nc'
			infile_name_anl = subprocess.check_output('ls '+analyses_dir+FILE_ANL, shell=True).strip()
			print ''
			print ' Analysis file:'
			print(infile_name_anl)
			print ''

			infile_ANL = netCDF4.Dataset(infile_name_anl,'r')

		     # Gather reanalysis vars
			rawLAT = infile_ANL.variables['g0_lat_0'][:]
			rawLON = infile_ANL.variables['g0_lon_1'][:]
			rawIVTU = infile_ANL.variables['VIEWVF_GDS0_EATM']
			rawIVTV = infile_ANL.variables['VINWVF_GDS0_EATM']
			rawIWV = infile_ANL.variables['TCWV_GDS0_SFC']

			idx_lat = np.abs( rawLAT - (lat_bot) ).argmin()  # find the nearest value
			if ( lon_bot >= 0): # its positive
				idx_lon = np.abs( rawLON - lon_bot ).argmin()
			else: # its negative
				idx_lon = np.abs( rawLON - (360 + lon_bot) ).argmin()

			theLAT  = rawLAT[idx_lat]
			theLON  = rawLON[idx_lon]
			print ' Closest lat/lon in analysis file:'
			print(str(theLAT)+'N & '+str(theLON)+'E')
			print ''

		     # Manipulation of the reanalysis data
			theIWV = rawIWV[idx_lat,idx_lon]
			theIVTU = rawIVTU[idx_lat,idx_lon]
			theIVTV = rawIVTV[idx_lat,idx_lon]
			theIVT = np.sqrt( theIVTU**2 + theIVTV**2 )
			print('IVT = '+str(theIVT))
			print('IWV = '+str(theIWV))
			print ''

		     # Write to list file
			#print('here-------------------------------------------')
			#line2write = occ_id +', '+ type_dataset +', '+ str(occ_year) +', '+ str(occ_month) +', '+ str(occ_day) + \
			#print(lat_round)
			#print(type(lat_round))
			line2write = occ_id +', '+ str(type_dataset) +', '+ occ_year +', '+ occ_month +', '+ occ_day + \
				', ' + str(occ_hour) +', '+ str(occ_minute) +', '+ str(occ_second) + \
				', ' + str(lat_bot) + ', ' + str(lon_bot) +  ', ' + str(hgt_obs_km_bot) + \
				', ' + str(year_round) +  ', ' + str(month_round) +  ', ' + str(day_round) +  ', ' + str(hour_round) + \
				', ' + str(lat_round) +  ', ' + str(lon_round) +  \
				', ' + str(theIVT) +  ', ' + str(theIWV) + '\n'
			#print(line2write)
			list_file.write(line2write)
			print('written to file')
			print('')


			#list_file.close()
			#sys.exit()

list_file.close()
print('')
print(' ====================')
print('   Graceful Stop')
print(' ====================')
print('')

