import numpy as np
import copy
import math


"""

Ian Ashpole, November 2022

-This is a series of functions used for creating the L3L and L3W datasets
that are described and analysed in Ashpole and Wiacek (2022, AMT, 	
https://doi.org/10.5194/amt-2022-90). The functions are:
    
    1) read_and_screen_L2_selected_variables
    2) get_L3L_L3W_from_L2
    3) save2nc_selected_variables

-each function contains a docstring and is heavily commented 

                                                        
"""



##########################################
# ----------------------------------------
# 1. read_and_screen_L2_selected_variables
# ----------------------------------------
##########################################



def read_and_screen_L2_selected_variables(l2file):

    """

    Purpose:
    -Read L2 data and screen following filters applied in L3 creation
    (as detailed in data users guide:
     https://www2.acom.ucar.edu/sites/default/files/documents/v8_users_guide_201812.pdf)
    -Returns screened L2 fields as a dictionary

    Input:
    l2file - an open hdf (.he5) MOP02 file

    Output:
    l2_data - a dictionary of screened L2 data, plus the lons and lats 
    that correspoind to each retrieval

    """


    # 1. get the following from l2 file:
    # ---a) latitudes and longitudes;
    # ---b) variables for quality checking data fields (ie screening);
    # ---c) data fields
    # (if variable names are not known, can use the below snippet)
    # for item in l2file['HDFEOS/SWATHS/MOP02/Data Fields'].items():
    #     print item
    l2_base_address='HDFEOS/SWATHS/MOP02/Data Fields/'

    # ---a)latitudes and longitudes
    l2_longitudes=l2file['/HDFEOS/SWATHS/MOP02/Geolocation Fields/Longitude'][:]
    l2_latitudes=l2file['/HDFEOS/SWATHS/MOP02/Geolocation Fields/Latitude'][:]

    # ---b) variables for performing quality checks on data fields (ie screening);
    sza=l2file[l2_base_address+'SolarZenithAngle'][:]
    pixel=l2file[l2_base_address+'SwathIndex'][:,0]
    l1re=l2file[l2_base_address+'Level1RadiancesandErrors'][:]

    # -----get radiance and radiance error values to calculate SNR
    # -----(from data users guide...)
    l1r=l1re[:,:,0]
    l1e=l1re[:,:,1]

    # -----get SNR for elements 5A and 6A, as these are where SNR is important
    # -----(used for screening later)
    l1_snr_5a=l1r[:,3]/l1e[:,3]
    l1_snr_6a=l1r[:,9]/l1e[:,9]

    # ---c) data fields
    l2_sfci=l2file[l2_base_address+'SurfaceIndex'][:].astype(float) #0=water, 1=land, 2=mixed
    l2_sfc_pressure=l2file[l2_base_address+'SurfacePressure'][:]
    l2_dry_air_column=l2file[l2_base_address+'DryAirColumn'][:]
    
    l2_tco_ret=l2file[l2_base_address+'RetrievedCOTotalColumn'][:,0]
    l2_tco_apr=l2file[l2_base_address+'APrioriCOTotalColumn'][:]
    l2_dfs=l2file[l2_base_address+'DegreesofFreedomforSignal'][:]
    l2_vmr_sfc_ret=l2file[l2_base_address+'RetrievedCOSurfaceMixingRatio'][:,0]
    l2_vmr_sfc_apr=l2file[l2_base_address+'APrioriCOSurfaceMixingRatio'][:,0]
    l2_ak_sfc_rowsum=l2file[l2_base_address+'AveragingKernelRowSums'][:,0]
    l2_ak_rowsums=l2file[l2_base_address+'AveragingKernelRowSums'][:,:]
  
    l2_ak_matrix=l2file[l2_base_address+'RetrievalAveragingKernelMatrix'][:]
    l2_ak_sfc_diagonal=l2_ak_matrix[:,0,0]
  
    
    # 2. calculate XCO (column averaged volume mixing ratio) 
    # from total column CO and dry air column
    tco_values=l2_tco_ret.astype('float64')
    dac_values=l2_dry_air_column.astype('float64')
    l2_xco=(tco_values/dac_values)*(10**9)
    
    
    # 3. screen the L2 data according to criteria used in L3 data creation,
    # as outlined in the data user's guide:

    # ---a) find invalid pixels 
    l2_invalid=np.where((sza >= 80) | (pixel == 3) | 
                        ((l1_snr_5a < 1000) & (l1_snr_6a < 400)))

    # ---b) set invalid pixels in l2 data fields to nan
    l2_sfci[l2_invalid]=np.nan
    l2_sfc_pressure[l2_invalid]=np.nan
    l2_tco_ret[l2_invalid]=np.nan
    l2_tco_apr[l2_invalid]=np.nan
    l2_dfs[l2_invalid]=np.nan
    l2_xco[l2_invalid]=np.nan
    l2_dry_air_column[l2_invalid]=np.nan
    l2_vmr_sfc_ret[l2_invalid]=np.nan
    l2_vmr_sfc_apr[l2_invalid]=np.nan
    l2_ak_sfc_rowsum[l2_invalid]=np.nan
    l2_ak_sfc_diagonal[l2_invalid]=np.nan
    

    # 4. record l2 data to dictionary 

    # ---a. create the dictionary
    l2_data={'l2_longitudes':l2_longitudes,
                    'l2_latitudes':l2_latitudes,
                    'l2_sfci':l2_sfci,
                    'sfc_pressure':l2_sfc_pressure,
                    'tco_ret':l2_tco_ret,
                    'tco_apr':l2_tco_apr,
                    'dfs':l2_dfs,
                    'xco':l2_xco,
                    'dry_air_column':l2_dry_air_column,
                    'vmr_sfc_ret':l2_vmr_sfc_ret,
                    'vmr_sfc_apr':l2_vmr_sfc_apr,
                    'ak_sfc_rowsum':l2_ak_sfc_rowsum,
                    'ak_sfc_diagonal':l2_ak_sfc_diagonal,
                    }
    
    
    # ---b. before outputting the l2_data, need to set any missing values
    # ---(currently flagged as -9999) to np.nan.
    for key in l2_data.keys():
        l2_data[key][l2_data[key] < -900]=np.nan
        
     
    # 5. NOTE that I found a small discrepancy in the outputted data fields
    # that it is necessary to correct for.
    # *********THE DISCREPANCY:**********
    # vmr_sfc values are given for every l2 retrieval, however, the 
    # retrieval "height" that the surface actually corresponds to varies 
    # due to topography. This is an issue because the retrieval 
    # has 10 levels: 9 are distributed evenly, every 100hPa between
    # 900 and 100hPa; and there is a 10th "surface" level, which is of 
    # variable thickness since the surface pressure varies. 
    # Where the surface elevation is above 900hPa (i.e. surface pressure is
    # lower), the corresponding layer becomes the surface level 
    # (e.g. if surface pressure is 825hPa, then the 900-800hPa level 
    # becomes the surface level (which = the 9th level)).
    # Correspondingly, the AK matrix for the 10th level is set to missing, 
    # and the surface level AK becomes the 9th level of the AK matrix, thus
    # creating a problem for the steps of the code so far. Essentially this
    # means that for l2 retrievals with a surface pressure < 900hPa, the 
    # AK rowsum and diagonal values for the surface level retrievals are 
    # currently missing.  
    # *********THE SOLUTION:**********
    # filling the AK info for surface level retrievals where surface 
    # pressure is < 900 hPa using info from the relevant level of the AK
    # matrix.
    
    # ----step 1: duplicate the fields for ak_sfc_rowsum and 
    # ak_sfc_diagonal: these will be "filled" with ak info from the new 
    # "surface" level in these instances       
    ak_sfc_rowsum_filled=copy.deepcopy(l2_ak_sfc_rowsum)
    ak_sfc_diagonal_filled=copy.deepcopy(l2_ak_sfc_diagonal)


    # ----step 2: identify instances whereby there is a value for 
    # vmr_sfc_ret but no corresponding ak_sfc_diagonal value: 
    # this indicates that the "surface" level is at a pressure below 900 hPa 
    # and therefore that the "surface" level in the ak matrix (the 10th 
    # retrievak level) is not the actual surface as the surface here is 
    # higher than that 
    # (the 10th retrieval level corresponds to pressures greater than 900hPa
    # and is "floating").
    # We call these "discrepancies"
    discrepancies=np.where(
        (np.isnan(l2_ak_sfc_diagonal)) & (~np.isnan(l2_vmr_sfc_ret))
        )
    
    
    # ----step 3: loop over the pressure levels to find the relevant 
    # retrieval level that is the new "surface"
    bottom_pressures=[900,800,700,600,500,400,300,200,100]
    level_indexes=[1,2,3,4,5,6,7,8,9] #for indexing the ak
    for bottom_pressure,level_index in zip(bottom_pressures,level_indexes):
    
        loi=np.where(
            (l2_sfc_pressure[discrepancies] < bottom_pressure) &
            (l2_sfc_pressure[discrepancies] >= bottom_pressure-100)
            )
            
        # ----step 4: if there are "surface" level retrievals at this 
        # pressure level, then get the corresponding AK rowsum and diagonal
        # value fields, screen them for data quality (using index defined
        # earlier), and copy the AK info to the "filled" AK surface fields
        if np.shape(loi)[1] > 0:

            ak_rowsum_thislevel=l2_ak_rowsums[:,level_index]
            ak_rowsum_thislevel[l2_invalid]=np.nan
            ak_sfc_rowsum_filled[discrepancies[0][loi]]=ak_rowsum_thislevel[discrepancies[0][loi]]
            
            ak_diagonal_thislevel=l2_ak_matrix[:,level_index,level_index]
            ak_diagonal_thislevel[l2_invalid]=np.nan
            ak_sfc_diagonal_filled[discrepancies[0][loi]]=ak_diagonal_thislevel[discrepancies[0][loi]]
        
            # end if 
            
        #end for (bottom_pressure,level index)

    
    # ----step 5: append the "filled" ak fields to the l2_data dictionary
    l2_data['ak_sfc_rowsum_filled']=ak_sfc_rowsum_filled
    l2_data['ak_sfc_diagonal_filled']=ak_sfc_diagonal_filled
        
    
    # 6. return (output) the l2 data dictionary
    return l2_data


    #finished




########################
# ----------------------
# 2. get_L3L_L3W_from_L2
# ----------------------
########################



def get_L3L_L3W_from_L2(l3file,l2_data):

    """

    Purpose:
    -Create 1x1deg L3L (land-only) and L3W (water-only) data from
    0.25deg L2 fields, using 1x1deg L3 data from the corresponding
    MOPITT time slot as the base (this follows the method outlined in
    Ashpole & Wiacek 2020 AMT and 2022 AMTD).
    -Returns dictionaries of the L3L and L3W data for this time slot.

    Input:
    l3file - an open hdf (.he5) MOP03 file for a given time slot
    l2_data - a dictionary of L2 retrievals for the corresponding time slot

    Output:
    l3l_data,l3w_data - dictionaries of 1x1 deg gridded l3l and l3w fields,
    with corresponding lons and lats from the l3 file, for this time slot
    
        
    """


    # 1. Specify L3 base address and read lat,lon & tco variables from
    # l3file. lat & lons will be used as iterables for L2 file, while 
    # tco will be used to determine which gridboxes to get data for 
    # (i.e. the non-nan gridboxes).
    # (to see data field shape: print file[address_to_variable].shape)
    l3_base_address='HDFEOS/GRIDS/MOP03/Data Fields/'
    l3_latitude=l3file[l3_base_address+'/Latitude/'][:]
    l3_longitude=l3file[l3_base_address+'/Longitude/'][:]
    l3_tco_ret=l3file[l3_base_address+'/RetrievedCOTotalColumnDay/'][:]


    # 2. Set missing elements in l3_tco_ret as np.nan and take a deepcopy of
    # this fld to use as the base for recording the surface-averaged L2 
    # data for each L3 gridbox to (l3l = l3_land_only, l3w = l3_water_only).
    l3_tco_ret[l3_tco_ret < -900]=np.nan

    # 3. initialise the l3l and l3w data dictionaries which are eventually
    # to be output by this function
    # (NOTE this step could probably be done using l2_data.keys(), but the 
    # current code is more obvious)
    l3l_data={'n_obs':copy.deepcopy(l3_tco_ret),
                    'tco_ret':copy.deepcopy(l3_tco_ret),
                    'tco_apr':copy.deepcopy(l3_tco_ret),
                    'dfs':copy.deepcopy(l3_tco_ret),
                    'xco':copy.deepcopy(l3_tco_ret),
                    'dry_air_column':copy.deepcopy(l3_tco_ret),
                    'sfc_pressure':copy.deepcopy(l3_tco_ret),
                    'vmr_sfc_ret':copy.deepcopy(l3_tco_ret),
                    'vmr_sfc_apr':copy.deepcopy(l3_tco_ret),
                    'ak_sfc_rowsum':copy.deepcopy(l3_tco_ret),
                    'ak_sfc_diagonal':copy.deepcopy(l3_tco_ret),
                    'ak_sfc_rowsum_filled':copy.deepcopy(l3_tco_ret),
                    'ak_sfc_diagonal_filled':copy.deepcopy(l3_tco_ret)
                    }
    
    l3w_data=copy.deepcopy(l3l_data)


    # 4. Get the index of non_nan (i.e. valid) l3 grid boxes, which is what
    # we create L3L and L3W for
    l3_non_nan=np.where(l3_tco_ret > -900.)

    
    # 5. specify a list of the vmr keys that appear as keys in the data 
    # dictionaries, because these need to be processed differently 
    # (the l2 retrievals that are used to create l3 are logged for vmrs,
    # but not for other variables)
    vmr_keys=['vmr_sfc_ret','vmr_sfc_apr','vmr_900_ret','vmr_900_apr',
              'vmr_800_ret','vmr_800_apr','vmr_600_ret','vmr_600_apr',
              'vmr_300_ret','vmr_300_apr']
     

    # 6. loop over l3_non_nan gridboxes and use the indices to get
    # longitude and latitude values for these gridboxes.
    for l3_gridbox in range(0,len(l3_non_nan[1])): #for 1 (l3_gridbox)


        # 7. Get corner coordinates for this l3 gridbox.
        # Remember, the l3 gridbox coordinates accessed from the
        # longitude,latitude vectors are the gridbox MID-POINT, therefore
        # need to get the lower and upper limits for these.
        l3lonoi_lower=math.floor(l3_longitude[l3_non_nan[0][l3_gridbox]])
        l3lonoi_upper=math.ceil(l3_longitude[l3_non_nan[0][l3_gridbox]])
        l3latoi_lower=math.floor(l3_latitude[l3_non_nan[1][l3_gridbox]])
        l3latoi_upper=math.ceil(l3_latitude[l3_non_nan[1][l3_gridbox]])


        # 8. Now find L2 pixels that fall within this 1x1 degree bounding box
        l2roi=np.where((l2_data['l2_longitudes'] >= l3lonoi_lower) &
            (l2_data['l2_longitudes'] <= l3lonoi_upper) &
            (l2_data['l2_latitudes'] >= l3latoi_lower) &
            (l2_data['l2_latitudes'] <= l3latoi_upper))


        # 9. get the L2 pixels inside l2roi that are over water (l2sfci = 0)
        # and land (l2sfci = 1)
        waterpixels=np.where(l2_data['l2_sfci'][l2roi] == 0)
        landpixels=np.where(l2_data['l2_sfci'][l2roi] == 1)


        # 10. record land and water ('L3L') and ('L3W') mean values 
        # from the bounded l2 data for this 1deg x 1deg l3 gridbox.
        # Use the l3l_data dictionary key (i.e. variable name) as the index,
        # as these are identical in l3 and l2 data dictionaries. This
        # greatly cuts down on lines of code!
        # NOTE: if the variable is not a vmr, then the average is simply the 
        # average of all bounded l2 data for this grid box. If it is a vmr 
        # however, the bounded data first need to be logged before they are 
        # averaged, and then the average needs to be antilogged to give the 
        # mean value to be recorded.
        for key in l3l_data.keys(): #for 2 (key)
            if key == 'n_obs':
                l3l_data[key][l3_non_nan[0][l3_gridbox]][l3_non_nan[1][l3_gridbox]]=np.sum(~np.isnan(l2_data['tco_ret'][l2roi][landpixels]))
                l3w_data[key][l3_non_nan[0][l3_gridbox]][l3_non_nan[1][l3_gridbox]]=np.sum(~np.isnan(l2_data['tco_ret'][l2roi][waterpixels]))

            elif key not in vmr_keys:
                l3l_data[key][l3_non_nan[0][l3_gridbox]][l3_non_nan[1][l3_gridbox]]=np.nanmean(l2_data[key][l2roi][landpixels])
                l3w_data[key][l3_non_nan[0][l3_gridbox]][l3_non_nan[1][l3_gridbox]]=np.nanmean(l2_data[key][l2roi][waterpixels])
            else:
                l3l_data[key][l3_non_nan[0][l3_gridbox]][l3_non_nan[1][l3_gridbox]]=np.exp(np.nanmean(np.log(l2_data[key][l2roi][landpixels])))
                l3w_data[key][l3_non_nan[0][l3_gridbox]][l3_non_nan[1][l3_gridbox]]=np.exp(np.nanmean(np.log(l2_data[key][l2roi][waterpixels])))
                
                #endelse
            #end for 2 (key)

        #end for 1 (l3_gridbox)


    # 11. Add l3 longitude and latitude fields to l3l_data and l3w_data
    # dictionaries
    l3l_data['longitude']=l3_longitude
    l3l_data['latitude']=l3_latitude
    l3w_data['longitude']=l3_longitude
    l3w_data['latitude']=l3_latitude
    

    # 12. Output the l3l_data and l3w_data dictionaries
    #print('completed get_L3L_L3W_from_L2')
    return l3l_data,l3w_data


    #finished




###############################
# -----------------------------
# 3. save2nc_selected_variables
# -----------------------------
###############################



def save2nc_selected_variables(data,savename):

    
    """
   
    Purpose:
    -creates (and saves!) a netcdf file from an input data dictionary
    -NOTE that the exact variables are required in the input data 
    dictionary, otherwise this function will crash. I have created 
    more flexible versions of the function, but this version is being
    used here for greatest clarity

    Input:
    data (type=dictionary). L3(-style) MOPITT data to be saved
    savename (type=string). Name (and location) to save the .nc file to.

    Output:
    netcdf file
    

    """

    
    # 0. import the netCDF4 library for working with netcdf
    import netCDF4 as nc4


    # 1. initialise the .nc file and create dimensions
    dataset=nc4.Dataset(savename,'w',format='NETCDF4_CLASSIC')

    latitude=dataset.createDimension('latitude',len(data['latitude']))
    longitude=dataset.createDimension('longitude',len(data['longitude']))
    time=dataset.createDimension('time', 1)


    # 2. create coordinate variables for file
    times=dataset.createVariable('time',np.float64, ('time',))
    latitudes=dataset.createVariable('latitude',np.float32,('latitude',))
    longitudes=dataset.createVariable('longitude',np.float32,('longitude',))


    # 3. create the actual data variables that we want to save
    n_obs=dataset.createVariable('n_obs',np.float32,('time','longitude','latitude'))
    dfs=dataset.createVariable('dfs',np.float32,('time','longitude','latitude'))
    tco_ret=dataset.createVariable('tco_ret',np.float32,('time','longitude','latitude'))
    tco_apr=dataset.createVariable('tco_apr',np.float32,('time','longitude','latitude'))
    xco=dataset.createVariable('xco',np.float32,('time','longitude','latitude'))
    dry_air_column=dataset.createVariable('dry_air_column',np.float32,('time','longitude','latitude'))
    sfc_pressure=dataset.createVariable('sfc_pressure',np.float32,('time','longitude','latitude'))
    vmr_sfc_ret=dataset.createVariable('vmr_sfc_ret',np.float32,('time','longitude','latitude'))
    vmr_sfc_apr=dataset.createVariable('vmr_sfc_apr',np.float32,('time','longitude','latitude'))
    ak_sfc_diagonal=dataset.createVariable('ak_sfc_diagonal',np.float32,('time','longitude','latitude'))
    ak_sfc_diagonal_filled=dataset.createVariable('ak_sfc_diagonal_filled',np.float32,('time','longitude','latitude'))
    ak_sfc_rowsum=dataset.createVariable('ak_sfc_rowsum',np.float32,('time','longitude','latitude'))
    ak_sfc_rowsum_filled=dataset.createVariable('ak_sfc_rowsum_filled',np.float32,('time','longitude','latitude'))

    #'sfci' is a special case that is present in original L3 data but not in
    #L3L and L3W that are created from L2 data. Test to see whether the key
    #'sfci' is present in the dictionary; only create the variable if it is.
    if 'sfci' in data:
        sfci=dataset.createVariable('sfci',np.float32,('time','longitude','latitude'))


    # 4. create the metadata
    dataset.description='MOPITT V8 L3-style data based on L2 NIR & TIR combined retrievals (MOP02J)'
    dataset.history='Created November 2022'
    dataset.source='Parent data downloaded from earthdata.nasa.gov, processed by Ian Ashpole @ SMU'
    dataset.fill_value=('-9999.')


    # 5. add variable attributes
    latitudes.units='degrees_north'
    longitudes.units='degrees_east'
    n_obs.units='count'
    dfs.units='unitless'
    tco_ret.units='mol/m2'
    tco_apr.units='mol/m2'
    xco.units='ppbv'
    dry_air_column.units='mol/m2'
    sfc_pressure.units='hPa'
    vmr_sfc_ret.units='ppbv'
    vmr_sfc_apr.units='ppbv'
    ak_sfc_diagonal.units='unitless'
    ak_sfc_diagonal_filled.units='unitless'
    ak_sfc_rowsum.units='unitless'
    ak_sfc_rowsum_filled.units='unitless'

    #'sfci' special case
    if 'sfci' in data:
        sfci.units='unitless'

    times.units='hours since 0001-01-01 00:00:00'
    times.calendar='gregorian'


    # 6. fill with data
    times[:]=1.
    latitudes[:]=data['latitude']
    longitudes[:]=data['longitude']
    n_obs[0,:,:]=data['n_obs']
    dfs[0,:,:]=data['dfs']
    tco_ret[0,:,:]=data['tco_ret']
    tco_apr[0,:,:]=data['tco_apr']
    xco[0,:,:]=data['xco']
    dry_air_column[0,:,:]=data['dry_air_column']
    sfc_pressure[0,:,:]=data['sfc_pressure']
    vmr_sfc_ret[0,:,:]=data['vmr_sfc_ret']
    vmr_sfc_apr[0,:,:]=data['vmr_sfc_apr']
    ak_sfc_diagonal[0,:,:]=data['ak_sfc_diagonal']
    ak_sfc_diagonal_filled[0,:,:]=data['ak_sfc_diagonal_filled']
    ak_sfc_rowsum[0,:,:]=data['ak_sfc_rowsum']
    ak_sfc_rowsum_filled[0,:,:]=data['ak_sfc_rowsum_filled']

    #'sfci' special case
    if 'sfci' in data:
        sfci[0,:,:]=data['sfci']



    # 7. close (write) the file
    dataset.close()


    #finished


 

#END
