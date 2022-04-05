from Used_modules import *
# Input file and different data set used in this project
#       Files :
#           - Cube : Velocity information from IRAM 30m data or any other data set
#           - Table : Polarization data from HAWC+/SOFIA NASA
#           - Table_super : Pol. data of HAWC+/SOFIA with superpixel method
#           - Herchel_map : Herchel space data set
#           - Planck_map  : Planck polarization data
#
#
#
#
#
'''Inputs: '''

cub = './data/C18O_combined_cube.fits'
table = './data/taurus_hawc.txt'
table = './data/L1495_superpixel.csv'
tablen = './data/L1495_superpixel_oct20.csv'

herschel_map = "./data/HGBS_tauN3_hires_column_density_map.fits"
planck_map = "./data/polang_ETaurus.fits"

region1 = [[64.5776885, 27.6264983], [64.5222854, 27.6271158],
           [64.5216194, 27.5747001], [64.5770037, 27.5724186]]
region2 = [[64.5311152, 27.5531077], [64.4952007, 27.6003805],
           [64.4559148, 27.5760812], [64.4974710, 27.5280108]]
region4 = [[64.5062131, 27.4997813], [64.4508000, 27.4878722],
           [64.4804781, 27.3790218], [64.5394474, 27.3948612]]
region3 = [[64.6620083, 27.4354328], [64.5696225, 27.5457689],
           [64.5144577, 27.5089470], [64.6078148, 27.3961557]]
mask_planck = [[66.0760873, 27.1380245], [64.5998007, 28.2190558], [
    64.2093623, 28.1877267], [63.5458751, 27.3486629], [65.3295163, 26.1248061]]
region3 = region = [[64.7006577, 27.4545927], [64.5813169, 27.5751658],
                    [64.4668327, 27.4831221], [64.5901029, 27.3592069]]
# region2=region=[[64.5487552,27.5581333],[64.5095419,27.6029979],[64.4491636,27.5457787],[64.4924545,27.5103427]]
# region1=region=[[64.6091896,27.6080715],[64.5708003,27.6529531],[64.5152572,27.6073431],[64.5528364,27.5646481]]

# regmimir=[[66.7250802,27.6739741],[65.1869945,28.7752978],[63.7636490,27.4210635],[65.4355444,26.3572498]]#first try ; je l ai donnee a Ph. pour verifier
regmimir = [[66.3642222, 28.2134334], [65.0625270, 28.5843339],
            [63.9593374, 27.2800004], [66.0369026, 26.2909323]]

n_H2 = 2*10**4  # volume density in this region region 1,2,3
density_cr = 21.825
plot_moment = 1
plot_moment_mimir = 1
plot_moment0 = 1
plot_stat_B = 0
plot_all_maps = 1
plot_spectra = 1
plot_stat = 1
plot_velocity = 1
plot_planck = 1
plot_her = 0
'''--------------------------'''

herschel = fits.open(get_pkg_data_filename(herschel_map))[0]
f = fits.getdata(herschel_map, ext=0)
c = f
f = np.log10(f)
(ys, xs) = f.shape


'''This to read and to plot the velocity disperssion  of inputed data'''

'''1. extract moment 0, 1, 2 from data'''
