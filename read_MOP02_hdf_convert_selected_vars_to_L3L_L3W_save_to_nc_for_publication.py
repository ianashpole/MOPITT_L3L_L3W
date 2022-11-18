import os
from os import system
from sys import exit
from datetime import datetime
from glob2 import glob as glob
import h5py
import resource,gc


"""
Ian Ashpole, Nov 2022

    
This code does the following:
-Open a list of MOP03 (MOPITT L3, as downloaded: 'L3O') files in hdf5 format
-Open the corresponding (by date) MOP02 (MOPITT L2) hdf file and read
    selected variables
-Using the L3 data and its 1x1deg grid as a base, create L3L (land-only)
    and L3W (water-only) data from the L2 fields
-Save daily L3L and L3W datasets as netcdf files
    missing values are set to 'nan' in the netcdf files as a postprocessing
    step using cdo


REQUIRES:
-function_create_L3L_and_L3W_from_MOP02_and_MOP03_hdf.py
-list of input MOP03 .he5 files
-list of input MOP02 .he5 files
-cdo installation (optional, for postprocessing)


USER REQUIRED EDITS:
-must specify the location of MOPITT L2 (MOP02) and L3 (MOP03) hdf files
    and also the location to save the output netcdf files to
    (locations are specified in Step 2)
    
    
"""


# (0. optional step for getting performance info)
starttime = datetime.now()
startusage=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


# 1. Import the custom functions for:
#--- a. processing MOP02 data
#--- b. creating L3L and L3W 
#--- c. saving the results to netcdf
from function_create_L3L_and_L3W_from_MOP02_and_MOP03_hdf import\
    read_and_screen_L2_selected_variables,\
        get_L3L_L3W_from_L2,\
            save2nc_selected_variables    


# 2. Specify data and save directories 
# (MOP02 is being processed; MOP03 used as base for creating L3L & L3W)
MOP03_datadir='path_to_L3_input_data'
MOP02_datadir='path_to_L2_input_data'
savedir='path_to_output_location'


# 3. Loop over MOP03 files & get file dt (using original L3 files ('L3O') 
# as base as only creating L3L & L3W for days where there is a L3O file)
flist=sorted(glob(MOP03_datadir+'MOP03J-20*.he5'))
nfiles=len(flist)
n_post_crash=0

for f in range(0,nfiles): # for 1 (f)
    dt=os.path.basename(flist[f]).split('-')[1]
    
    if int(dt) >= 20010901: # if 1 (dt, put in place in case of code crash to avoid reprocessing the same files)
        print('file',f,'of',nfiles,';',dt)
        n_post_crash=n_post_crash+1
    
    
        # 4. Open MOP03 hdf5.
        l3file=h5py.File(flist[f],'r')
    
    
        # 5. Search for MOP02 file with the same dt as the L3 file and
        # open that too. If this L2 file it doesn't exist, there is an 
        # error that needs resolving.
        l2_flist=glob(MOP02_datadir+'MOP02J-'+dt+'*.he5')
    
        if len(l2_flist) == 0:
            print('ERROR: no L2 file corresponding to L3 file for ',dt)
            exit()
        else:
            l2file=h5py.File(l2_flist[0])
    
    
        # 6. Pass l2file to function that reads selected variables and screens
        # them according to L3 data creation guidelines outlined in data user
        # guides. This returns a dictionary of screened L2 data.
        l2_data=read_and_screen_L2_selected_variables(l2file)
        
        
        # 7. pass the l3file and l2_data dictionary to function that creates
        # l3l and l3w from the screened l2_data using the l3file as a base.
        l3l_data,l3w_data=get_L3L_L3W_from_L2(l3file,l2_data)
        
        
        # 8. Close the l2 and l3 files to clear memory & speed up code    
        l2file.close()
        l3file.close()
        
        
        # 9. Pass the data dictionaries to function to save the data as 
        # netcdf files. 
        # NOTE 1: Need to pass this function a savename.
        # NOTE 2: after the .nc file is created, postprocess to set 
        # missing values to nan. cdo is the easiest tool for doing this
        # (this step can just be commented out if cdo is not installed
        # or is not seen as necessary; NOTE that if this is done, 
        # 'savename' needs to be given to the save2nc function call in 
        # place of 'foo.nc')
        savename=savedir+'L3L/MOPITT_v8.L3L.from_MOPO2J.selected_variables.'+dt+'.nc'
        save2nc_selected_variables(l3l_data,savedir+'foo.nc')
        system('cdo -s setmissval,nan '+savedir+'foo.nc '+savename)
    
        savename=savedir+'L3W/MOPITT_v8.L3W.from_MOPO2J.selected_variables.'+dt+'.nc'
        save2nc_selected_variables(l3w_data,savedir+'foo.nc')
        system('cdo -s setmissval,nan '+savedir+'foo.nc '+savename)
        
        
        # 10. (optional)
        # Delete arrays that were created during loop and clear them from 
        # system memory to speed up code
        del(l2_data,l3l_data,l3w_data,l2file,l3file)
        gc.collect()
        
        
        print('mean time per file:',(datetime.now()-starttime)/float(n_post_crash))
        
        
            
        #end if 1 (dt, put in place in case of code crash to avoid reprocessing the same files))
        
 
    #end for 1 (f)


#finished (optional performance metrics printout below)
print(startusage,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
print('Finished @ (hh.mm.ss.ms): {}'.format(datetime.now()))
print('Time elapsed (hh.mm.ss.ms): {}'.format(datetime.now()-starttime))
exit()
