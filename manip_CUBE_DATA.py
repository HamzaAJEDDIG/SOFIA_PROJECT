from Used_modules import *
def extract_moments(cub):
    '''this will return moment 0,1,2 in selected region
        '''
    hi_data = fits.open(cub)  # Open the FITS file for reading
    cube = SpectralCube.read(hi_data)  # Initiate a SpectralCube
    hi_data.close()
    # sub_cub=cube[51:106,:,:] # This to select just where there is detection
    sub_cub = cube.spectral_slab(4.5 * u.km / u.s, 7 * u.km / u.s)
    moment_0 = sub_cub.with_spectral_unit(u.km/u.s).moment(order=0)  # Zero-th moment
    moment_1 = sub_cub.with_spectral_unit(u.km/u.s).moment(order=1)  # the first moment
    moment_2 = sub_cub.with_spectral_unit(u.km/u.s).moment(order=2)  # the seconde moment
    linewidth_sigma = sub_cub.with_spectral_unit(u.km/u.s).linewidth_sigma()
    return cube, sub_cub, moment_0, moment_1, moment_2, linewidth_sigma
