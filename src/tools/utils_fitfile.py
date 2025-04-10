from astropy.io import fits

def readfit(path):
    
    with fits.open(path) as hdu:
        # https://docs.astropy.org/en/stable/io/fits/
        dict_fitfile = {
            'data'       : hdu[0].data, # .astype(np.uint8)
            'freqs'      : hdu[1].data  ['Frequency'][0], # extract frequency axis in MHz 
            'timeax'     : hdu[1].data  ['Time'][0],      # extract time axis (in seconds, with timeax[0] = 0)
            'dT'         : hdu[0].header['CDELT1'],       # extract time resolution (in seconds)
            'date'       : hdu[0].header['DATE-OBS'],     # take first file
            'T0'         : hdu[0].header['TIME-OBS'],     # take first file
            'instrument' : hdu[0].header['INSTRUME'],
            'content'    : hdu[0].header['CONTENT']}
        
    return dict_fitfile