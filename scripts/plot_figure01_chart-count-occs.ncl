;________________________________________________________________
;
;   NCL script to plot a map using the Robinson rojection
;    and overlay the location of the lowest tangent point
;    for each observed GNSS RO profile.
;    The RO obs are read from raw IODA (JEDI) netCDF4 files.
;
;   Reference for Robinson labeling:
;    https://docs.dkrz.de/doc/visualization/sw/ncl/examples/source_code/dkrz-ncl-add-latitude-longitude-annotation-to-robinson-projection-example.html
;
;    ~ mjm 2023-10-20
;________________________________________________________________

; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   

begin

function add_robinson_labels(wks,map,latspc,lonspc)
local lat_values, lon_values, lat_labels, lon_labels, \
      nlat, nlon, txres, x_in, y_in, x_out, y_out
begin
  minlat =  -90
  maxlat =   90
  minlon = -180
  maxlon =  180

;-- pick some "nice" values for the latitude labels
  lat_values = ispan(minlat,maxlat,latspc) * 1.
  lon_values = ispan(minlon,maxlon,lonspc) * 1.
  nlat       = dimsizes(lat_values)
  nlon       = dimsizes(lon_values)

;-- create the labels; add space before right labels, and after left labels
  lat_labels = where(lat_values.lt.0,abs(lat_values)+"~S~o~N~S",lat_values+"~S~o~N~N")
  lat_labels = where(lat_values.eq.0,"0~S~o~N~",lat_labels)
  lat_labels = where(abs(lat_values).eq.90,"",lat_labels) ;-- no label at lat=abs(90)

;---------------------
;-- latitude labels
;---------------------
  txres                  =  True     ;-- set text resources
  txres@txFontHeightF    =  0.010
  txres@txFontThicknessF =  2

  y_in  = lat_values                         ;-- plot coordinates
  x_in  = new(dimsizes(y_in),float)          ;-- plot coordinates
  x_out = new(dimsizes(x_in),float)          ;-- for NDC coordinates
  y_out = new(dimsizes(y_in),float)          ;-- for NDC coordinates

;-- left latitude labels
  txres@txJust = "CenterRight"
  x_in = 0.0
  ;x_in = -180.0
  datatondc(map,x_in,y_in,x_out,y_out)
  gsn_text_ndc(wks, lat_labels, x_out-0.012, y_out, txres)

;-- right latitude labels
  txres@txJust = "CenterLeft"
  x_in = -1.0
  ;x_in = 180.0
  datatondc(map,x_in,y_in,x_out,y_out)
  gsn_text_ndc(wks, lat_labels, x_out+0.012, y_out, txres)

;---------------------
;-- longitude labels
;---------------------
  lon_labels = where(lon_values.lt.0,abs(lon_values)+"~S~o~N~W",lon_values+"~S~o~N~E")
  lon_labels = where(lon_values.eq.0,"",lon_labels) ;-- no label at lon=abs(0)
  ;lon_labels = where(lon_values.eq.0,"0~S~o~N~",lon_labels)
  lon_labels = where(abs(lon_values).eq.180,"180~S~o~N~",lon_labels) 
  ;lon_labels = where(abs(lon_values).eq.180,"",lon_labels) ;-- no label at lon=abs(180)

;-- longitude labels
  txres@txJust           = "CenterCenter"
  txres@txFontHeightF    =  0.009
  txres@txFontThicknessF =  2

  x_in  := lon_values
  y_in  := new(dimsizes(x_in),float)
  x_out := new(dimsizes(x_in),float)
  y_out := new(dimsizes(y_in),float)

;-- add top longitude labels
  y_in = 90.0
  datatondc(map,x_in,y_in,x_out,y_out)
  gsn_text_ndc(wks, lon_labels, x_out, y_out+0.012, txres)

;-- add bottom longitude labels
  y_in = -90.0
  datatondc(map,x_in,y_in,x_out,y_out)
  gsn_text_ndc(wks, lon_labels, x_out, y_out-0.012, txres)

;-- that's it
  return(map)

end   ;-- end function add_robinson_labels



;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Begin main script
;;;;;;;;;;;;;;;;;;;;;;;;;;


 ; User input
	;root_path = "../"
	root_path = "/discover/nobackup/mjmurph4/"
	;date_str = "2021-12-04"
	;date_str = "2021-12-12"
	date_str = "2022-02-15"
	;time_str = "1800"
	lat_int = 30 ; lat label interval
	lon_int = 30 ; lon lable interval

	constel_sel = "common-ro" ; i.e. non-COSMIC2 & non-Spire
	;constel_sel = "only-cosmic2"
	use_spire = False ; use the file containing CSDA Spire RO?
	;constel_sel = "only-spire"
	;use_spire = True ; use the file containing CSDA Spire RO?

	;type = "x11"
	type = "pdf"

	;threshold_occ_bottom_km = 100 ; in km
	;threshold_occ_bottom_km = 30 ; in km
	;threshold_occ_bottom_km = 20 ; in km
	;threshold_occ_bottom_km = 10 ; in km
	threshold_occ_bottom_km = 5 ; in km
	;threshold_occ_bottom_km = 3 ; in km
	;threshold_occ_time_offset = 6 ; in hours, zero is right at analysis time


	output_path = "../plots/"
	;output_path = root_path + "plots/"
	if ( use_spire ) then
		infile_name_root = "gdas1_spnasa"
		;outfile_name_suffix = "only-spire"
	else
		infile_name_root = "gdas1"
		;outfile_name_suffix = "no-spire"
	end if


 ; Time manipulation
	campaign_year_str = str_get_cols(date_str,0,3)
	campaign_mon_str = str_get_cols(date_str,5,6)
	campaign_day_str = str_get_cols(date_str,8,9)

	date_tmp = str_sub_str(date_str,"-","")
	date_str_fmt = date_str
	date_filename_ioda = date_tmp + "_????z"


 ; Read in the IODA file
	path_ioda_root = root_path + \
			  "projects/hiaper/mjmurphy/modeling/input_forcing/ioda_format/gdas_bufr/gnssro/" + \
			  "2022/??/??/"
			  ;"processed_data/"
;/discover/nobackup/mjmurph4/projects/hiaper/mjmurphy/modeling/input_forcing/gdas_bufr/2022/02/15/

	PATH = path_ioda_root 
	;PATH = path_ioda_root + campaign_year_str +"/"+ campaign_mon_str +"/"+ campaign_day_str +"/"
	FILE = infile_name_root+"."+date_filename_ioda+".gpsro.tm00_ioda.nc4"
	FILE = "ioda_ncep_gnssro_pt6h_"+date_filename_ioda+".nc4"
	;FILE = "ioda_ncep_gnssro_pt6h_20220215_??00z.nc4"
	infile_str = systemfunc("ls "+PATH+FILE)

	print("")
	print("IODA obs file to read: " + infile_str)
	print("")

	infiles = addfiles(infile_str,"r")
	;ListSetType(infiles, "join") ; doesn't work on JEDI output files
	;ListSetType(infiles, "cat") ; default


 ; Gather the data we need from IODA file
	rawLAT  = infiles[:]->latitude
	rawLON  = infiles[:]->longitude
	rawRECN = infiles[:]->sequenceNumber ; differentiates occs
	rawTYPE = infiles[:]->satelliteIdentifier ; RO sat ID
	rawHGT  = infiles[:]->height ; geometric height in meters
	;rawTime = infiles[:]->time ; time offset to analysis time, in hours
	;printVarSummary(rawLAT)


     ; Read the dims from each file
	; KLUDGE to work aournd addfiles not working with JEDI output
	nfiles = dimsizes(infile_str)
	nlocs_file = new(nfiles,integer)
	do ifile=0,nfiles-1
		infile = addfile(infile_str(ifile),"r")
		nlocs_file(ifile) = getfiledimsizes(infile)
		delete(infile)
	end do
	;print(nlocs_file)


 ; Manipultion of the IODA data
	nlocs = dimsizes(rawRECN)

     ; Stagger the record_numbers for each file
	; KLUDGE to work aournd addfiles not working with JEDI output
	idx_end_stag = new(nfiles,integer)
	do ifile=0,nfiles-1
		if ( ifile .eq. 0) then
			idx_end_stag(ifile) = 0 +  nlocs_file(ifile) -1
		else
			idx_end_stag(ifile) = idx_end_stag(ifile-1) + nlocs_file(ifile)
		end if
	end do
	;print(idx_end_stag)

	rawRECN_stag = new(nlocs,typeof(rawRECN))
	offset_stag = 100000 ; has to more than nlocs
	do ifile=0,nfiles-1
		if ( ifile .eq. 0 ) then
			idx_start = 0
		else
			idx_start = idx_end_stag(ifile-1)
		end if
		
		idx_end = idx_end_stag(ifile)
		;print(idx_start+" & " +idx_end)

		rawRECN_stag(idx_start:idx_end) = rawRECN(idx_start:idx_end) + offset_stag
		offset_stag = offset_stag + offset_stag
	end do
	;print(rawRECN +" ==> "+ rawRECN_stag)


     ; Get unique occ identifier
	unique_id_occ = get_unique_values(rawRECN_stag) ; record number seems to do the trick
	printVarSummary(unique_id_occ)
	;print(unique_id_occ)

	;count_missing_occ = num(ismissing(unique_id_occ))
	;print(count_missing_occ)

	noccs = dimsizes(unique_id_occ)
	occ_low_lat  = new(noccs,float)
	occ_low_lon  = new(noccs,float)
	occ_low_alt  = new(noccs,float)
	occ_low_type = new(noccs,float)

	do iocc = 0,noccs-1
		idx_occ_tmp = ind( rawRECN_stag .eq. unique_id_occ(iocc))

		;subsetTime = rawTime(idx_occ_tmp)
		subsetTYPE = rawTYPE(idx_occ_tmp)
		subsetLAT  = rawLAT(idx_occ_tmp)
		subsetLON  = rawLON(idx_occ_tmp)
		subsetHGT  = rawHGT(idx_occ_tmp) /1000. ; convert meters to km 

		if ( subsetHGT(0) .gt. threshold_occ_bottom_km ) then
		;if ( (subsetHGT(0) .gt. threshold_occ_bottom_km) .or. \
		;     (abs(subsetTime(0)) .gt. threshold_occ_time_offset) ) then
			delete( [/ idx_occ_tmp /] )
			delete( [/ subsetLAT, subsetLON, subsetHGT, subsetTYPE /] )
			;delete( [/ subsetLAT, subsetLON, subsetHGT, subsetTime, subsetTYPE /] )

			continue ; skip it

		else
			;print(subsetHGT(0))
			;print(subsetTime(0))
			occ_low_type(iocc) = subsetTYPE(0) ; same for all heights
			occ_low_lat(iocc)  = subsetLAT(0) ; first is lowest		
			occ_low_lon(iocc)  = subsetLON(0) ; first is lowest		
			occ_low_alt(iocc)  = subsetHGT(0) ; first is lowest

			delete( [/ idx_occ_tmp /] )
			delete( [/ subsetLAT, subsetLON, subsetHGT, subsetTYPE /] )
			;delete( [/ subsetLAT, subsetLON, subsetHGT, subsetTime, subsetTYPE /] )
			;delete( [//] )
		end if		


	end do

	;print(occ_low_lat)
	;print(occ_low_lon)
	;print(occ_low_alt)
	;print(occ_low_type)


     ; Subset by constellation
	outfile_name_suffix = constel_sel

	if ( constel_sel .eq. "common-ro")
		idx_selected = ind( ( (occ_low_type .lt. 750) .or. (occ_low_type .gt. 755) ) \
					.and. (occ_low_type .ne. 269) \
					.and. (occ_low_type .ne. 3) \
					.and. (occ_low_type .ne. 4) \
					.and. (occ_low_type .ne. 5) ) ; not COSMIC-2 & not Spire & not MetOp
		;idx_selected = ind( ( (occ_low_type .lt. 750) .or. (occ_low_type .gt. 755) ) \
		;			.and. (occ_low_type .ne. 269) ) ; not COSMIC-2 & not Spire
		idx_selected_02 = ind( (occ_low_type .ge. 3) .or. (occ_low_type .le. 5) ) ; MetOp
		;occ_color = "forestgreen"
		occ_color = "blue"
		occ_symbol_02 = 5 ; Xs

	else if ( constel_sel .eq. "only-cosmic2")
		idx_selected = ind( (occ_low_type .ge. 750) .and. (occ_low_type .le. 755) ) ; only COSMIC-2
		occ_color = "purple"

	else if ( constel_sel .eq. "only-spire")
		idx_selected = ind( (occ_low_type .eq. 269) )
		occ_color = "blue"

	else
		print("")
		print(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
		print("  Constellation not recognized")
		print("   therefore ABORTING now!!")
		print(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
		print("")
		exit

	end if
	end if
	end if

	occ_low_lat_selected = occ_low_lat(idx_selected)
	occ_low_lon_selected = occ_low_lon(idx_selected)
	;print(occ_low_type(idx_selected))
	occs_selected = dimsizes(occ_low_lat_selected)
	;print(occ_low_type(idx_selected))

     ; Confirm the right constellations were selected
	constel_id_used =  get_unique_values(occ_low_type(idx_selected))
	print("")
	print("Constellation IDs included:")
	print( constel_id_used )
	print("")
	;exit

     ; Counting
	count_missing_occ = num(ismissing(occ_low_alt))
	count_used_occ = noccs - count_missing_occ
	print("")
	print("  Total Occulations = "+tostring(noccs))
	print("")
	print("  Threshold for lowest tp = "+threshold_occ_bottom_km +" km")
	print("  Occulations within thresholds = "+tostring(count_used_occ))
	print("  Selected RO within thresholds  = "+tostring(occs_selected))
	;print("  Non-Cosmic2 within thresholds = "+tostring(occs_noncos2))
	;print("  Occulations outside thresholds = "+tostring(count_missing_occ))
	print("")



 ; We generate plots, but what kind do we prefer?
	oufile_name = output_path + "chart_occ-count-lowest_"+date_str+"_thres"+threshold_occ_bottom_km+"km_"+ \
				outfile_name_suffix+"."+type
	wks = gsn_open_wks(type,oufile_name)


 ; Setting up the map
	mpres = True
	;mpres@gsnMaximize   		= True      	; Maximize plot in frame

	mpres@gsnDraw               =  False
	mpres@gsnFrame              =  False
	mpres@gsnLeftString         = ""
	mpres@gsnRightString        = ""

	mpres@vpWidthF              =  0.9
	mpres@vpHeightF             =  0.8
	mpres@vpXF                  =  0.05
	mpres@vpYF                  =  0.95

	;mpres@cnFillOn              =  True
	;mpres@cnLineLabelsOn        =  False

	mpres@mpPerimOn 		=  False
	mpres@mpProjection 		= "Robinson"
	mpres@mpCenterLonF              = 180
	;mpres@mpCenterLonF              = 160


     ; Grid Line Options
	mpres@mpGridAndLimbOn        	= True      	; turn on lat/lon lines
	mpres@mpGridLineColor        	= "Black" 	; set grid color
	;mpres@mpGridLineColor        	= "Gray" 	; set grid color
	;mpres@mpGridLineColor        	= "grey70" 	; set grid color
	mpres@mpGridLatSpacingF 	=  lat_int
	mpres@mpGridLonSpacingF 	=  lon_int
	mpres@mpGridLineDashPattern  	= 2          	; lat/lon lines as dashed
	;mpres@mpGridLineThicknessF  	= 0.7  		; Sets lat & lon gridline thickness
	;mpres@mpGridSpacingF         	= 15  		; Sets lat & lon gridline spacing
	;mpres@mpGridSpacingF         	= 2  		; Sets lat & lon gridline spacing
	;mpres@mpGridMaskMode 		= "MaskNotOcean"; only show grid lines over ocean
	;mpres@mpPerimOn             	= False     	; turn off perimeter


 ; Plotting
	plot = gsn_csm_map(wks,mpres)        ; draw global map


     ; Attach latitude labels using function
	map  = add_robinson_labels(wks,plot,lat_int,lon_int)


     ; Add RO lowest tangent point to map
	opts_poly1 = True
	;opts_poly1@gsMarkerIndex = 4
	;opts_poly1@gsMarkerSizeF	= 0.004 
	opts_poly1@gsMarkerIndex = 1
	;opts_poly1@gsMarkerThicknessF 	= 3    	; thickness of makers
	;opts_poly1@gsMarkerSizeF	= 0.006 
	opts_poly1@gsMarkerSizeF	= 0.011 
	opts_poly1@gsMarkerColor = occ_color
	poly_mark_selected = gsn_add_polymarker(wks,plot,occ_low_lon_selected,occ_low_lat_selected,opts_poly1)

	if ( constel_sel .eq. "common-ro")
		opts_poly2 = True
		opts_poly2@gsMarkerIndex = occ_symbol_02
		opts_poly2@gsMarkerSizeF = 0.005
		;opts_poly2@gsMarkerSizeF = 0.002
		opts_poly2@gsMarkerColor = occ_color
		occ_low_lat_selected_02 = occ_low_lat(idx_selected_02)
		occ_low_lon_selected_02 = occ_low_lon(idx_selected_02)
		poly_mark_selected_02 = gsn_add_polymarker(wks,plot,occ_low_lon_selected,occ_low_lat_selected_02,opts_poly2)
	end if

     ; Add desired region outlines to map
	opts_polyline_tropics = True
	opts_polyline_tropics@gsLineColor = "red"     ; color of lines
	opts_polyline_tropics@gsLineThicknessF = 2.5           ; thickness of lines
	;opts_polyline_tropics@gsLineThicknessF = 1.5           ; thickness of lines

	lat_tropics_max = 20.0
	lat_tropics_min = -20.0

	lat_tropics_max_line = (/lat_tropics_max, lat_tropics_max, lat_tropics_max/)
	lat_tropics_min_line = (/lat_tropics_min, lat_tropics_min, lat_tropics_min/)
	lon_tropics_line = (/ 0, 180, 360/)

	poly_tropics_top = gsn_add_polyline(wks,plot,lon_tropics_line,lat_tropics_max_line,opts_polyline_tropics)
	poly_tropics_bot = gsn_add_polyline(wks,plot,lon_tropics_line,lat_tropics_min_line,opts_polyline_tropics)


	opts_polyline_nepac = True
	opts_polyline_nepac@gsLineColor = "forestgreen"     ; color of lines
	opts_polyline_nepac@gsLineThicknessF = 2.5           ; thickness of lines
	;opts_polyline_nepac@gsLineThicknessF = 1.5           ; thickness of lines

	lat_nepac_max = 55.0
	lat_nepac_min = 20.0
	lon_nepac_max = 360 - 125.0
	lon_nepac_min = 360 - 179.0

	lat_nepac_max_line = (/lat_nepac_max, lat_nepac_max, lat_nepac_max/)
	lat_nepac_min_line = (/lat_nepac_min, lat_nepac_min, lat_nepac_min/)
	lon_nepac_line = (/ lon_nepac_min, lon_nepac_max/)

	lat_nepac_line = (/ lat_nepac_min, lat_nepac_max/)
	lon_nepac_max_line = (/lon_nepac_max, lon_nepac_max, lon_nepac_max/)
	lon_nepac_min_line = (/lon_nepac_min, lon_nepac_min, lon_nepac_min/)

	poly_nepac_top = gsn_add_polyline(wks,plot,lon_nepac_line,lat_nepac_max_line,opts_polyline_nepac)
	poly_nepac_bot = gsn_add_polyline(wks,plot,lon_nepac_line,lat_nepac_min_line,opts_polyline_nepac)
	poly_nepac_left = gsn_add_polyline(wks,plot,lon_nepac_min_line,lat_nepac_line,opts_polyline_nepac)
	poly_nepac_right = gsn_add_polyline(wks,plot,lon_nepac_max_line,lat_nepac_line,opts_polyline_nepac)

	draw(plot)
	frame(wks)
	;maximize_output(wks,False)

end
exit 

