# MOPITT_L3L_L3W
Code for creating the MOPITT L3L and L3W datasets described in Ashpole and Wiacek 2022 (AMT: https://doi.org/10.5194/amt-2022-90)  
&nbsp;  
## Repository contents
This repository contains 2 python code files:
1. read_MOP02_hdf_convert_selected_vars_to_L3L_L3W_save_to_nc.py
    - a wrapper code for specifying in/out data locations, opening files and passing them into functions in #2 for processing. Also contains optional postprocessing calls to cdo.
2. function_create_L3L_and_L3W_from_MOP02_and_MOP03_hdf.py
    - 3 functions to 
        1. process the input L2 data; 
        2. create L3L and L3W; and 
        3. save them to netcdf
  
## Requirements
- Python packages required: numpy, h5py, netCDF4
- Optional additional packages required: cdo (a command-line tool used here for postprocessing the netcdf files that are created) 
    - (NOTE that the required cdo functionality can easily be replicated in python)  

## Source data locations
- MOPITT L3 files: https://doi.org/10.5067/TERRA/MOPITT/MOP03J_L3.008
- MOPITT L2 files: https://doi.org/10.5067/TERRA/MOPITT/MOP02J_L2.008 

Both addresses above direct to MOPITT Version 8 joint thermal- & near-infrared (TIR-NIR) data. The code in this repository will work with TIR-only and NIR-only files, as well as other data versions, although it has only been tested up to V9. Alternative MOPITT datasets can be searched for at https://search.earthdata.nasa.gov (L2/L3 MOPITT files have the tags "MOP02" and "MOP03", resp.)

## User edits required to run the code
- Data in/out locations need to be specified in `read_MOP02_hdf_convert_selected_vars_to_L3L_L3W_save_to_nc.py` 
    - Optionally, cdo functionality can be toggled on/off