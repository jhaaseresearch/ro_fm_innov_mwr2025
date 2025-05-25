#________________________________________________________________
#
#   Python script to read JEDI hofx output in its native
#    netCDF format and calculate statistics on RO obs. 
#    This example is set up for plotting the observation
#    minus the background (OmB).
#    
#   ~ mjm 2022-06-29
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

nargs_req = 6
if ( nargs == nargs_req ):
	arrecon_year_str = sys.argv[1]
	iop_str = sys.argv[2]
	series_id = sys.argv[3]
	exp_id = sys.argv[4]
	region = sys.argv[5]
	output_root = sys.argv[6]
else:
	print("")
	print("")
	print("   Wrong number of arguements, ABORTING!!")
	print("")
	print("   You gave " +str(nargs)+ " of the "+str(nargs_req)+" required args" )
	print("")
	print("     Usage:")
	print("       " + script_name + " <arrecon_year> <iop_str> <series_id> <exp_id> <region> <output_root>")
	print("")
	print("     Example:")
	print("       python " + script_name + " 2022 full 05 05 nepac ../")
	print("")
	#print("     NOTE: the stats arg can be")
	#print("            1) full, 2) inside-ar, or 3) outside-ar")
	print("")
	print("")

	sys.exit()



 # User input
input_root = '/discover/nobackup/mjmurph4'
input_root_path = input_root + '/projects/hiaper/mjmurphy/modeling/output_fv3jedi/'
#aoutput_root = '../'
#iop_num = int(iop_str[-2:]) # last 2 characters are iop  number
#threshold_top = 15 # in km, cuttoff for top of data
#threshold_top = 20 # in km, cuttoff for top of data
threshold_top = 30 # in km, cuttoff for top of data
inc_hgt_bins_m = 200
#inc_hgt_bins_m = 250
#inc_hgt_bins_m = 500


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


# -- strip info from filename
year_str = arrecon_year_str
#PATH_OBS = input_root_path + '/'+year_str+iop_str+'/series'+series_id+'/'+exp_id+'/'+year_str+'-01-11_????/hofx/'
#PATH_OBS = input_root_path + '/'+year_str+iop_str+'/series'+series_id+'/'+exp_id+'/'+year_str+'-02-??_????/hofx/' # only February
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
hgt_bins_upper_m = range(0+inc_hgt_bins_m,threshold_top_m+inc_hgt_bins_m,inc_hgt_bins_m)
#hgt_bins_upper_m = range(0+inc_hgt_bins_m,20000+inc_hgt_bins_m,inc_hgt_bins_m)
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
hgt_km_occ_list = []
ba_omb_list = []
ba_omb_percent_list = []
#ba_err_list = []
#ba_err_eff_list = []
#ba_err_percent_list = []
ba_omb_binned_list = []
ba_omb_percent_binned_list = []
#ba_err_binned_list = []
#ba_err_eff_binned_list = []
#ba_err_percent_binned_list = []


 # Loop through all the files to get occultations
for ifile in infile_name_obs_list:
	print('')
	print('=======================================================================================')

	infile_name_obs = ifile
	#print(infile_name_obs)
	infile_OBS = h5py.File(infile_name_obs,'r')
	#infile_OBS = netCDF4.Dataset(infile_name_obs,'r')
	print('')
	print(' Observation file:')
	print('   ' + infile_name_obs)


     # Gather the data we need
	ba_hofx_raw = infile_OBS['hofx']['bendingAngle'][:]
	#ba_hofx_raw = infile_group_hofx.variables['bending_angle'][:]
	#print(ba_hofx_raw)

	ba_obs_raw  = infile_OBS['ObsValue']['bendingAngle'][:]
	#ba_err_raw  = infile_OBS['ObsError']['bendingAngle'][:]
	#ba_err_eff_raw  = infile_OBS['EffectiveError']['bendingAngle'][:]

	hgt_obs_raw = infile_OBS['MetaData']['height'][:] # Geometric altitude in m
	lat_obs_raw = infile_OBS['MetaData']['latitude'][:] 
	lon_obs_raw = infile_OBS['MetaData']['longitude'][:] # west is negative
	#id_sat_raw  = infile_group_meta.variables['occulting_sat_id'][:] # west is negative
	#id_ref_raw  = infile_group_meta.variables['reference_sat_id'][:] # west is negative
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
	#id_ref_tmp  = id_ref_raw[index_hgt]

	ba_obs_tmp  = ba_obs_raw[index_hgt]
	ba_hofx_tmp = ba_hofx_raw[index_hgt]
	#ba_err_tmp  = ba_err_raw[index_hgt]
	#ba_err_eff_tmp  = ba_err_eff_raw[index_hgt]
	
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

     # Mask the missing values
	#ba_hofx_missing_value = -3.3687953e+38 # discovered by trial and error
	ba_hofx_mask = np.ma.masked_where( ba_hofx == ba_hofx_missing_value, ba_hofx )
	#ba_err_eff_mask = np.ma.masked_where( ba_err_eff == ba_hofx_missing_value, ba_hofx )
	#sys.exit()

     # Calculate the main variables
	ba_omb = (ba_obs - ba_hofx_mask ) * 1000 # convert to mrad
	ba_omb_percent = ( (ba_obs - ba_hofx_mask) / ba_hofx_mask ) * 100
	hgt_obs_km  = hgt_obs / 1000
	#ba_err_percent = ( (max(ba_err) - ba_err ) / max(ba_err) ) * 100
	#ba_err_eff = ba_err_eff_mask

	#nheights = ba_obs.shape
	#nheights = nheights[0]
	#print(nheights)

	#mean_omb = np.mean(ba_omb) # result is nan

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
			#print(lon_top)
			print(lon_bot)

			hgt_obs_km_sel = hgt_obs_km[rec_num == iid]
			ba_omb_sel = ba_omb[rec_num == iid]
			ba_omb_percent_sel = ba_omb_percent[rec_num == iid]
			#ba_err_sel = ba_err[rec_num == iid]
			#ba_err_eff_sel = ba_err_eff[rec_num == iid]
			#ba_err_percent_sel = ba_err_percent[rec_num == iid]

		      # Bin the data for the given occultation
			#   (taking mean of values in the given occ over the height range of the given bin)
			ba_omb_binned = np.ma.masked_all(shape=(nbins),dtype='float128') # initialize a fully masked array
			ba_omb_percent_binned = np.ma.masked_all(shape=(nbins),dtype='float128') # initialize a fully masked array
			#ba_omb_binned = np.ma.masked_all(shape=(nbins)) # initialize a fully masked array
			#ba_omb_percent_binned = np.ma.masked_all(shape=(nbins)) # initialize a fully masked array
			#ba_err_binned = np.ma.masked_all(shape=(nbins)) # initialize a fully masked array
			#ba_err_eff_binned = np.ma.masked_all(shape=(nbins)) # initialize a fully masked array
			#ba_err_percent_binned = np.ma.masked_all(shape=(nbins)) # initialize a fully masked array
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
				#ba_err_tmp = ba_err_sel[idx_hgt_bin_occ]
				#ba_err_eff_tmp = ba_err_eff_sel[idx_hgt_bin_occ]
				#ba_err_percent_tmp = ba_err_percent_sel[idx_hgt_bin_occ]
		
				ndims_tmp = ba_omb_percent_tmp.shape[0]
				if (  ndims_tmp == 0 ): # just stop here if its empty
					del(ba_omb_percent_tmp)
					continue

				ba_omb_binned[ibin] = np.ma.mean(ba_omb_tmp, axis=0)
				ba_omb_percent_binned[ibin] = np.ma.mean(ba_omb_percent_tmp, axis=0)
				#ba_err_binned[ibin] = np.ma.mean(ba_err_tmp, axis=0)
				#ba_err_eff_binned[ibin] = np.ma.mean(ba_err_eff_tmp, axis=0)
				#ba_err_percent_binned[ibin] = np.ma.mean(ba_err_percent_tmp, axis=0)

				del(ba_omb_tmp)
				del(ba_omb_percent_tmp)
				#del(ba_err_tmp)
				#del(ba_err_percent_tmp)


		      # Save the important vars in lists for each occ
			hgt_km_occ_list.append( hgt_obs_km[rec_num == iid] )
			ba_omb_list.append( ba_omb[rec_num == iid] )
			ba_omb_percent_list.append( ba_omb_percent[rec_num == iid] )
			#ba_err_list.append( ba_err[rec_num == iid] )
			#ba_err_eff_list.append( ba_err_eff[rec_num == iid] )
			#ba_err_percent_list.append( ba_err_percent[rec_num == iid] )

			ba_omb_binned_list.append(ba_omb_binned)
			ba_omb_percent_binned_list.append(ba_omb_percent_binned)
			#ba_err_binned_list.append(ba_err_binned)
			#ba_err_eff_binned_list.append(ba_err_eff_binned)
			#ba_err_percent_binned_list.append(ba_err_percent_binned)


noccs = len(hgt_km_occ_list)
print('')
print('================================')
print('   Total number of occs = ' + str(noccs))
#print(hgt_km_occ_list)


 # -- Count values & Square of the OmB values for calculation of RMSE later
 # --  combining this because its probably more effecient
ba_omb_binned_squared_list = []
ba_omb_percent_binned_squared_list = []
count_missing_bins = np.ma.zeros(shape=(nbins)) # initialize an array full of zeros
for iocc in range(noccs):
     # First square the OmB values for each occ
	ba_omb_binned_squared_list.append(np.ma.power(ba_omb_binned_list[iocc],2))
	ba_omb_percent_binned_squared_list.append(np.ma.power(ba_omb_percent_binned_list[iocc],2))

     # Next count the missing values for each occ and each bin
	for ibin in range(nbins):
		count_missing_bins[ibin] = count_missing_bins[ibin] + \
						np.ma.count_masked(ba_omb_binned_list[iocc][ibin], axis=None)

count_obs_bins = noccs - count_missing_bins
#print(refrac_count_missing_bins)
print("Number of occs per bin:")
print(count_obs_bins)
print('================================')
print('')


 # Square of the OmB values for calculation of RMSE later
#ba_omb_binned_squared_list = []
#ba_omb_percent_binned_squared_list = []
#for iocc in range(noccs): # must do this for each occ in list or else the np.ma part fails
#	ba_omb_binned_squared_list.append(np.ma.power(ba_omb_binned_list[iocc],2))
#	ba_omb_percent_binned_squared_list.append(np.ma.power(ba_omb_percent_binned_list[iocc],2))


 # -- Calculate the mean stats
ba_omb_mean = np.ma.mean(ba_omb_binned_list, axis=0)
ba_omb_std = np.ma.std(ba_omb_binned_list, axis=0)
ba_omb_rmse = np.ma.sqrt(np.ma.mean(ba_omb_binned_squared_list, axis=0))
#ba_omb_rmse = np.ma.sqrt(np.ma.mean(ba_omb_binned_list, axis=0)**2)
ba_omb_percent_mean = np.ma.mean(ba_omb_percent_binned_list, axis=0)
ba_omb_percent_std = np.ma.std(ba_omb_percent_binned_list, axis=0)
ba_omb_percent_rmse = np.ma.sqrt(np.ma.mean(ba_omb_percent_binned_squared_list, axis=0))
#ba_omb_percent_rmse = np.ma.sqrt(np.ma.mean(ba_omb_percent_binned_list, axis=0)**2)



 # Write out to ascii
#outfile_header = "Ellipsoid_ht_bin_mid(km) mean_omb(mrad) mean_omb_percent(%) stddev_omb(mrad) stddev_omb_percent(%) nobs(#/bin)"
outfile_header = "altitude_bin_mid(km) mean_omb(mrad) mean_omb_percent(%) stddev_omb(mrad) stddev_omb_percent(%) rmse_omb(mrad) rmse_omb_percent(%) nobs(#/bin)"
#outfile_header = "Ellipsoid_ht_bin_mid(km) mean_refrac_diff(N-units) mean_refrac_diff_percent(%) stddev_refrac_diff(N-units) stddev_refrac_diff_percent(%) nobs(#/bin)"
#output_array = zip(hgt_upper_km_out, refrac_diff_mean, refrac_diff_percent_mean, \
output_array = zip(hgt_mid_km, ba_omb_mean, ba_omb_percent_mean, \
                        ba_omb_std, ba_omb_percent_std, \
                        ba_omb_rmse, ba_omb_percent_rmse, \
                        count_obs_bins)
print(outfile_header)
print(output_array)
print(hgt_mid_km)
print(ba_omb_mean)
print(ba_omb_percent_mean)
print(ba_omb_std)
print(ba_omb_percent_std)
print(count_obs_bins)

output_path_ascii = output_root + 'stats_ascii/'
if not os.path.exists(output_path_ascii):
	os.makedirs(output_path_ascii)

outfile_name = output_path_ascii + \
		'stats_fv3jedi-omb_arrecon'+year_str+iop_str+ \
			'_exp'+series_id+exp_id+'_bin-'+str(inc_hgt_bins_m)+'m_'+region+'.txt'

np.savetxt(outfile_name, list(output_array), delimiter=" ", header=outfile_header, fmt='%1.4f')

print('')
print(' Stats file (ascii) created:')
print('   ' + outfile_name)
print('')



sys.exit()


