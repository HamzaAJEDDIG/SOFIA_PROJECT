import matplotlib.pyplot as plt
from Used_modules import *
from Input import *
from manip_CUBE_DATA import *
print(" ---------------------- main project SOFIA------------------------")
print(" -                                                               -")
print(" -                        H. Ajeddig                             -")
print(" -                       CEA - Saclay                            -")
print(" -                                                               -")
print(" -                                                               -")
print(" -                                                               -")
print(" -                        Have fun !                             -")
# Manipulation CUBE data
cube, sub_cub, moment_0, moment_1, moment_2, linewidth_sigma = extract_moments(cub)
print(moment_0.shape)
#
'''1./ First lets manipulate, reproject, transfrom, clean the data  
        '''
from manip_SOFIA_DATA import *
x, y, ra, dec, pa_s, sig_pa_s, deg_s, sig_deg_s, vlr, vlr_disp, Sig_linew = sofia_dat(
    table, moment_1, moment_2, linewidth_sigma)
#print( "x coordinates", x)
#print( "y coordinates", y)

''' 2./ Multi-Guassian fits of spectra of IRAM 30m data 
    '''
from plot_spectra import *
plot_spectra(cube, sub_cub, moment_0, moment_1, moment_2, linewidth_sigma,x, y, ra, dec, pa_s, sig_pa_s, deg_s, sig_deg_s, vlr, vlr_disp, Sig_linew)
#plt.show()

'''3./ Over plot different data sets and match their coordinates
    '''
from plot_overlap import *
plot_moment0()

'''4./ Calculate moments : see H. Ajeddig PhD manuscript 
    '''
from plot_moments import *
plot_moments()

'''Here, I added a new data set from another telescope for comparaison and claculating the error
    '''
from plot_momets_mimir import *
plot_moment_mimir()

'''5./ At this stage, we will reproject all the data to one plot and derive some statistics 
    '''
from reproject_data import *
px_v, py_v, Vlr, vlr_d, sigma_linew, pa_so, sig_pa_so, deg_so, sig_deg_so = map_data(
    herschel_map, ra, dec, pa_s, sig_pa_s, deg_s, sig_deg_s, vlr, vlr_disp, Sig_linew)
'''6./ Show the data at large scale and put vectors in the plot
    '''
from plot_all_data import *
px_v, py_v, Vlr, vlr_d, sigma_linew, pa_so, sig_pa_so, deg_so, sig_deg_so,p_a =plot_all_maps(px_v, py_v, Vlr, vlr_d, sigma_linew, pa_so, sig_pa_so, deg_so, sig_deg_so)

'''7./ The most difficult part, calculating the DCF method using all the data sets
    '''
from DCF_claculations import *
VLR3, dis_VLR3=plot_velocity()
print(len(VLR3))
B, V, planck_Bf, hawc = DFC_sub_region(region, Vlr, vlr_d, sigma_linew, pa_so, sig_pa_so, deg_so, sig_deg_so, VLR3, dis_VLR3,p_a,VLR3, dis_VLR3)

'''8./ Comparing the results with literature, show histogram and DCF results 
    '''
from plot_Planck_data import *
plot_planck(p_a,px_v, py_v)
plt.show()