"""
Description:

    Check which data has not been process yet and return the filenames that need
    to be processed.

    @author: Andrea F. Battaglia (andrea.francesco.battaglia@irsol.usi.ch)

    2024/12/12: first version
"""
import os

def return_fits2process(directory_fitfiles, directory_output):

   ### list all folders in dir_proc
   folders_proc = os.listdir(directory_output)
   folders_proc.sort()   

   ### List all raw FITS files in dir_raw that start with "meteoswiss" and ends with "fit"
   files_raw = [f for f in os.listdir(directory_fitfiles) if f.startswith('meteoswiss') and f.endswith('fit')]
   files_raw.sort()

   ### From flies_raw, extract the date of the observation
   dates_raw = [f.split('_')[1] for f in files_raw]

   ### Find the indices of the files in dates_raw that have no corresponding folder in folders_proc
   idx = [i for i, d in enumerate(dates_raw) if d not in folders_proc]

   ### Get the files that have no corresponding folder in folders_proc
   files2process = [files_raw[i] for i in idx]
   return files2process
