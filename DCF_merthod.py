import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import astropy.units as u
from astropy.io import fits
from spectral_cube import SpectralCube
from astroquery.esasky import ESASky
from astroquery.utils import TableList
from astropy.wcs import WCS
from reproject import reproject_interp
from numpy import *
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from collections import OrderedDict
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.data import get_pkg_data_filename
import aplpy
from scipy.signal import find_peaks

from scipy.stats import iqr
###############################################################################################################
#                               This to compute the magnetic field stregth using DCF method                   #
#                               the input :                                                                   #
#                                          1/. fits file of cub data containing velocity                      #
#                                          2/. fits file conataining the backgroud map                        #
#                                            (in this case Herscel data)                                      #
#                                          3/. table containing polarization data organize                    #
#                                                RA, DEC ,PA,  SIF_PA,  PDEG,  SIG_PDEG                       #
#                                          4/. fits file containing Planck data to be displayed (PA)          #
###############################################################################################################
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


cube, sub_cub, moment_0, moment_1, moment_2, linewidth_sigma = extract_moments(cub)

'''This to plot cub data and to see moment 0, 1, 2 and linewidth_sigma maps
        '''


def sofia_dat(table, moment_1, moment_2, linewidth_sigma):
    '''Input :This will read sofia data if they are in tabe format organized as
        ra, dec , I, sigI, P, sigP, d, sigd
        table must be in txt format to this stage to be allowed to use this function
        output : ra, dec in degrees and value of polarization angle and degrees with their incertitudes
        '''
    csv = np.genfromtxt(tablen, delimiter=",")
    # pa_s : is the B-field angle
    ra = csv[1:, 2]
    dec = csv[1:, 3]
    pa_s = csv[1:, 8]
    sig_pa_s = csv[1:, 9]
    deg_s = csv[1:, 6]
    sig_deg_s = csv[1:, 7]
    '''this is to extract data from cub where there is a polarization mesearement in sofia data '''
    w = WCS(moment_1.header)
    r, d = w.all_world2pix(np.array(ra), np.array(dec), 1)
    radec = []
    # here to estimate velocity component
    Vel = moment_1.hdu.data
    (ys, xs) = moment_1.hdu.data.shape
    var = moment_2.hdu.data
    sig_v = np.sqrt(moment_2.hdu.data)
    sig_linew = (linewidth_sigma.hdu.data)
    print( len(d))
    r = r.astype(int)
    d = d.astype(int)
    vlr = []
    vlr_disp = []
    Sig_linew = []
    x = []
    y = []
    for i in range(len(d)):
        radec.append([int(r[i]), int(d[i])])
    for j in range(len(radec)):
        #print radec[j][0],radec[j][1],Vel[radec[j][1],radec[j][0]],sig_v[radec[j][1],radec[j][0]]
        vlr.append(Vel[radec[j][1], radec[j][0]])
        vlr_disp.append(sig_v[radec[j][1], radec[j][0]])
        Sig_linew.append(sig_linew[radec[j][1], radec[j][0]])
        x.append(int(radec[j][0]))
        y.append(int(radec[j][1]))

    return x, y, ra, dec, pa_s, sig_pa_s, deg_s, sig_deg_s, vlr, vlr_disp, Sig_linew


x, y, ra, dec, pa_s, sig_pa_s, deg_s, sig_deg_s, vlr, vlr_disp, Sig_linew = sofia_dat(
    table, moment_1, moment_2, linewidth_sigma)
# check----------------------------------------------------------------------------------------
'''def sofia_dat(table,moment_1,moment_2,linewidth_sigma):
    Input :This will read sofia data if they are in tabe format organized as
        ra, dec , I, sigI, P, sigP, d, sigd
        table must be in txt format to this stage to be allowed to use this function
        output : ra, dec in degrees and value of polarization angle and degrees with their incertitudes

    csvn = np.genfromtxt('L1495_superpixel_oct20.csv', delimiter=",")
    #pa_s : is the B-field angle
    ran          =   csvn[1:,2]
    decn         =   csvn[1:,3]
    pa_sn        =   csvn[1:,8]
    sig_pa_sn    =   csvn[1:,9]
    deg_sn       =   csvn[1:,6]
    sig_deg_sn   =   csvn[1:,7]
    this is to extract data from cub where there is a polarization mesearement in sofia data
    wn = WCS(moment_1.header)
    rn,dn=wn.all_world2pix(np.array(ran),np.array(decn),1)
    radecn=[]
    #here to estimate velocity component
    Veln=moment_1.hdu.data
    (ysn,xsn)=moment_1.hdu.data.shape
    varn=moment_2.hdu.data
    sig_vn=np.sqrt(moment_2.hdu.data)
    sig_linewn=(linewidth_sigma.hdu.data)
    print len(dn)
    rn=rn.astype(int)
    dn=dn.astype(int)
    vlrn=[]
    vlr_dispn=[]
    Sig_linewn=[]
    xn=[]
    yn=[]
    for i in range(len(dn)):
        radecn.append([int(rn[i]),int(dn[i])])
    for j in range(len(radecn)):
        #print radec[j][0],radec[j][1],Vel[radec[j][1],radec[j][0]],sig_v[radec[j][1],radec[j][0]]
        vlrn.append(Veln[radecn[j][1],radecn[j][0]])
        vlr_dispn.append(sig_vn[radecn[j][1],radecn[j][0]])
        Sig_linewn.append(sig_linewn[radecn[j][1],radecn[j][0]])
        xn.append(int(radecn[j][0]))
        yn.append(int(radecn[j][1]))

    return xn,yn,ran,decn,pa_sn,sig_pa_sn,deg_sn,sig_deg_sn,vlrn, vlr_dispn, Sig_linewn

xn,yn,ran,decn,pa_sn,sig_pa_sn,deg_sn,sig_deg_sn,vlrn, vlr_dispn,Sig_linewn=sofia_dat(tablen,moment_1,moment_2,linewidth_sigma)
'''  # ----------------------------------------------------------------------------------------------
print( "x coordinates", x)
print( "y coordinates", y)
if plot_spectra == 1:
    # region=[[64.6620083,27.4354328],[64.5696225,27.5457689],[64.5144577,27.5089470],[64.6078148,27.3961557]]
    w = WCS(moment_0.header)
    radec = []
    for i in range(4):
        r, d = w.all_world2pix(region[i][0], region[i][1], 1)
        radec.append([int(r), int(d)])
        mask = np.zeros((ys, xs), dtype=np.float)
        polys = np.array([radec])
        import cv2
        import os

    for i in polys.tolist():
        cv2.fillConvexPoly(mask, np.array(i), 1)
        mask[mask == 0] = np.nan
    row, col = 3, 3
    loop = 1
    comp = 0
    num = 0
    fig3 = plt.figure(num, figsize=(14, 8))

    y_sp = []
    print( 'number of points vlr ', len(vlr))
    for i in range(len(vlr)):
        if mask[y[i], x[i]] == 1:
            if loop == row*col+1:
                num += 1
                fig3 = plt.figure(num, figsize=(14, 8))
                fig3 = plt.figure(figsize=(14, 8))
                fig3.patch.set_facecolor('white')
                loop = 1
                axes = fig3.add_subplot(row, col, loop)
            axes = fig3.add_subplot(row, col, loop)
            #plt.subplots_adjust(wspace=0, hspace=0.0, top=0.95, bottom=0.05, left=0.17, right=0.845)
            x_specter = cube.with_spectral_unit(u.km/u.s).spectral_axis
            #print x_specter.value
            np.exp(x_specter.value)
            #print type(x_specter)
            y_specter = cube[:, x[i], y[i]]

            peaks, _ = find_peaks(y_specter, height=1., width=np.std(x_specter.value))
            #print peaks
            plt.plot(x_specter, y_specter, drawstyle='steps-mid',
                     linewidth=2.5, solid_joinstyle='round')
            plt.plot(x_specter[peaks], y_specter[peaks], "o")
            #print x[i],y[i]
            y_sp.append(y_specter)
            #sub_cub[:, x[i],y[i]].quicklook()

            def gaus(x, a, x0, sigma):
                return a*exp(-(x-x0)**2/(2*sigma**2))

            #popt,pcov = curve_fit(gaus,x,y,p0=[1,mean,sigma])
            #print popt,pcov
            xpeak = x_specter[peaks].value
            ypeak = y_specter[peaks].value
            # vlr_pic.append(np.array(xpeak))
            for i in range(len(ypeak)):
                # vlr_pic.append(xpeak[i])
                # int_pic.append(ypeak[i])
                a = ypeak[i]
                x0 = xpeak[i]
                # sigma=np.std(xpeak)
                sigma = 0.2
                #print a,x0

                b2 = gaus(x_specter.value, a, x0, sigma)
                #print b2
                plt.plot(x_specter.value, b2, label='fit')
            comp += 1
            loop += 1


'''    plt.figure(10)
    print comp

    plt.plot(x_specter,np.average(y_sp),drawstyle='steps-mid',linewidth=2.5,solid_joinstyle='round')
    plt.xlabel('Mean of centeroid velocity (km/s)')
    plt.ylabel('Mean of Ta (K)')
    plt.title('Region 4')
    xl1=[4.5,4.5,4.5,4.5]
    yl1=[-0.2,1,1.5,2]
    plt.plot(xl1,yl1,'--',color='red')
    xl2=[7,7,7,7]
    yl2=[-0.2,1,1.5,2]
    plt.plot(xl2,yl2,'--',color='red')
    xl3=[2,4,6,10]
    yl3=[0,0,0,0]
    plt.ylim(np.min(y_sp/comp),0.5)
    plt.plot(xl3,yl3,'--',color='green',linewidth=1.5,alpha=0.8


    print "number of point ",comp'''
if plot_moment0 == 1:
    # region 3
    # 41 29 169.480922985 -16.2049814923
    # 59 57 169.376872037 -16.1858810959
    # 54 47 169.411658723 -16.1959396596
    # 56 49 169.400171648 -16.1962653938

    # region 1
    # 36 26 169.499331827 -16.2014251961
    # 44 28 169.476690096 -16.2131161666
    # 34 36 169.481247197 -16.1743289799
    # 57 51 169.395581799 -16.1933013677
    # 64 61 169.355205006 -16.1889098907
    # region 2
    # 72 63 169.332551848 -16.2005779879
    # 64 76 169.319739134 -16.1567580176
    # 82 76 169.280680553 -16.1964784861
    # 62 74 169.331233799 -16.1564437546

    #x_all=[21, 24, 21, 24, 24, 26, 26, 29, 31, 34, 26, 29, 31, 34, 36, 39, 41, 26, 29, 31, 36, 69, 71, 44, 76, 36, 39, 41, 44, 46, 76, 79, 36, 39, 76, 79, 39, 76, 79, 81, 39, 44, 46, 49, 76, 41, 44, 46, 49, 41, 49, 51, 76, 44, 49, 51, 79, 49, 51, 54, 46, 49, 51, 54, 56, 76, 46, 49, 51, 56, 59, 54, 56, 59, 61, 59, 61, 64, 49, 59, 61, 71, 69, 71, 74, 76, 71, 74, 76, 79, 71, 74, 74, 76, 56, 59, 61, 76, 56, 59, 61, 56, 59, 61, 81, 84, 86, 56, 59, 61, 81, 84, 61, 64, 81, 84, 64, 84, 86, 89]
    #y_all=[2, 2, 4, 4, 7, 9, 12, 12, 12, 12, 14, 14, 14, 14, 14, 14, 14, 17, 17, 17, 17, 22, 24, 27, 27, 29, 29, 29, 29, 29, 29, 29, 32, 32, 32, 32, 34, 34, 34, 34, 37, 37, 37, 37, 37, 39, 39, 39, 39, 42, 42, 42, 42, 44, 44, 44, 44, 47, 47, 47, 49, 49, 49, 49, 49, 49, 52, 52, 52, 52, 52, 54, 54, 54, 54, 57, 57, 57, 59, 59, 59, 69, 72, 72, 72, 72, 74, 74, 74, 74, 77, 77, 79, 79, 82, 82, 82, 82, 84, 84, 84, 87, 87, 87, 87, 87, 87, 89, 89, 89, 89, 89, 92, 92, 92, 92, 94, 94, 94, 97]
    xl = [52, 57, 54, 64, 36, 44, 34, 57, 64, 72, 64, 82, 62]
    yl = [84, 91, 94, 91, 26, 28, 36, 51, 61, 63, 76, 76, 74]

    hi_datafile = './data/C18O_combined_cube.fits'
    hi_data = fits.open(hi_datafile)  # Open the FITS file for reading
    cube = SpectralCube.read(hi_data)  # Initiate a SpectralCube
    hi_data.close()

    cub = cube.spectral_slab(4.5 * u.km / u.s, 7 * u.km / u.s)
    sub_cube_slab = cub
    print (cub.shape)
    _, b, _ = cube.world[0, :, 0]  # extract latitude world coordinates from cube
    _, _, l = cube.world[0, 0, :]  # extract longitude world coordinates from cube
    moment_0 = sub_cube_slab.with_spectral_unit(u.km/u.s).moment(order=0)
    print( "Moment 0")
    print( "max and min", np.nanmax(moment_0.hdu.data), np.nanmin(moment_0.hdu.data))
    print('Moment_0 has units of: ', moment_0.unit)
    cub_sn = cube.spectral_slab(0 * u.km / u.s, 2. * u.km / u.s)

    sn_moment_0 = cub_sn.with_spectral_unit(u.km/u.s).moment(order=0)
    SN = 3*np.nanstd(sn_moment_0.hdu.data)
    print (SN, np.nanstd(sn_moment_0.hdu.data), np.nanmean(sn_moment_0.hdu.data))
    moment_0.hdu.data[np.where(moment_0.hdu.data < SN)] = 0
    fig = plt.figure(figsize=(12, 8))
    row, col = 2, 2
    loop = 1

    ax = fig.add_subplot(111, projection=WCS(cub[0].header))
    im = ax.imshow(moment_0.hdu.data, cmap='rainbow')
    ax.invert_yaxis()  # Flips the Y axis

    ax.set_xlabel("RA (J2000)", fontsize=10)
    ax.set_ylabel("DEC (J2000)", fontsize=10)

    cont_levels1 = [np.nanmax(moment_0.hdu.data)*i for i in (0.3, 0.5, 0.7)]
    CS2 = ax.contour(moment_0.hdu.data, levels=cont_levels1, colors='black', linewidths=1)
    plt.scatter(x, y, s=40, c='white', marker='X', edgecolors='green')
    cax = plt.axes([0.78, 0.13, 0.03, 0.74])
    cbar = plt.colorbar(im, pad=.07, cax=cax)
    cbar.set_label('I (K. km.$s^{-1}$)', size=10)

    labels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13']
    #ax = fig.add_subplot(111, projection=moment_0.wcs)
    #plt.subplots_adjust(bottom = 0.1)
    plt.scatter(xl, yl, marker='o', c='black', s=100, edgecolors='white', alpha=0.4)

    for label, i, j in zip(labels, xl, yl):
        plt.annotate(
            label,
            xy=(i, j), xytext=(-20, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.4', fc='yellow', alpha=0.5),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    # plt.scatter([57,59],[94,94],color='blue')
    # plt.scatter([57,59],[51,51],color='blue')

if plot_moment == 1:
    # Initiate a figure and axis object with WCS projection information
    fig1 = plt.figure(figsize=(12, 8))
    row, col = 2, 2
    loop = 1
    sig_v = np.sqrt(moment_2)
    moment = [moment_0, moment_1, moment_2, linewidth_sigma]
    title = ['$C^{18}O$ Integrated intensity from 4.5 km.$s^{-1}$ to 7 km.$s^{-1}$',
             'Centeroid velocity map', 'Velocity variance map', '1D velocity disperssion ']

    for i in range(4):

        ax = fig1.add_subplot(row, col, loop, projection=moment[i].wcs)
        im = ax.imshow(moment[i].hdu.data, cmap='rainbow')
        cbar = plt.colorbar(im, pad=.07)
        if i == 0:
            im.set_clim(-0.5, 2)
            cbar.set_label('I (K. km.$s^{-1}$)', size=10)
        if i == 1:
            im.set_clim(4.5, 7)
            cbar.set_label('$V_{LSR}$ (K. km.$s^{-1}$)', size=10)
        if i == 2:
            im.set_clim(0.1, 0.5)
            cbar.set_label('$V_{LSR}$ ($km^2$.$s^{-2}$)', size=10)
        if i == 3:
            im.set_clim(0.2, 0.8)
            cbar.set_label('$\sigma _{V_{LSR}}$ ( km.$s^{-1}$)', size=10)

        plt.title(title[i])
        ax.invert_yaxis()  # Flips the Y axis

        ax.set_xlabel("RA (J2000)", fontsize=10)
        ax.set_ylabel("DEC (J2000)", fontsize=10)
        # Add a colorbar

        # Overplot column density contours
        cont_levels1 = [np.nanmax(moment[0].hdu.data)*i for i in (0.3, 0.5, 0.7)]
        CS2 = ax.contour(moment[0].hdu.data, levels=cont_levels1, colors='black', linewidths=1)
        plt.scatter(x, y, s=20, c='white', marker='X', edgecolors='green', alpha=1)
        #plt.scatter(xn,yn,s=30,c = 'black', marker = 'o',edgecolors = 'red',alpha=1)

        loop += 1
plot_moment_mimir = 1
if plot_moment_mimir == 1:
    # Initiate a figure and axis object with WCS projection information
    cubm = './data/taurus_filament_13co.fits'
    hi_datam = fits.open(cubm)  # Open the FITS file for reading
    cubem = SpectralCube.read(hi_datam)  # Initiate a SpectralCube
    hi_datam.close()
    # sub_cub=cube[51:106,:,:] # This to select just where there is detection
    sub_cubm = cubem.spectral_slab(4. * u.km / u.s, 9 * u.km / u.s)
    moment_0m = sub_cubm.with_spectral_unit(u.km/u.s).moment(order=0)  # Zero-th moment
    moment_1m = sub_cubm.with_spectral_unit(u.km/u.s).moment(order=1)  # the first moment
    moment_2m = sub_cubm.with_spectral_unit(u.km/u.s).moment(order=2)
    linewidth_sigmam = sub_cubm.with_spectral_unit(u.km/u.s).linewidth_sigma()
    # SpectralCube.linewidth_fwhm
    fig1m = plt.figure(figsize=(12, 8))
    rowm, colm = 2, 2
    loopm = 1
    sig_vm = np.sqrt(moment_2m)
    #region=[[64.6620083,27.4354328],[64.5696225,27.5457689],[64.5144577,27.5089470],[64.6078148,27.3961557]] #big region
    #regmimir=[[66.3736266,28.0036378],[65.0246842,28.7203045],[64.6501804,27.5980822],[66.3202391,26.6980821]] #small region
    (ysm, xsm) = moment_0m.hdu.data.shape
    wm = WCS(moment_0m.header)
    radecm = []
    for i in range(4):
        rm, dm = wm.all_world2pix(regmimir[i][0], regmimir[i][1], 1)
        radecm.append([int(rm), int(dm)])
        maskm = np.zeros((ysm, xsm), dtype=np.float)
        polysm = np.array([radecm])
        import cv2
        import os

    for im in polysm.tolist():
        cv2.fillConvexPoly(maskm, np.array(im), 3)
        maskm[maskm == 0] = np.nan

    # moment_1m.hdu.data[np.where(moment_0m.hdu.data<SN)]=0
    # moment_2m.hdu.data[np.where(moment_0m.hdu.data<SN)]=0
    # linewidth_sigmam.hdu.data[np.where(moment_0m.hdu.data<SN)]=0
    momentm = [moment_0m, moment_1m, moment_2m, linewidth_sigmam]
    title = ['$C^{18}O$ Integrated intensity from 4. km.$s^{-1}$ to 9 km.$s^{-1}$',
             'Centeroid velocity map', 'Velocity variance map', '1D velocity disperssion ']

    for i in range(4):

        axm = fig1m.add_subplot(rowm, colm, loopm, projection=momentm[i].wcs)
        imm = axm.imshow(momentm[i].hdu.data, cmap='rainbow')
        cbarm = plt.colorbar(imm, pad=.07)
        axm.imshow(maskm, cmap='gray',alpha=0.5)
        if i == 0:
            imm.set_clim(-0.5, 5)
            cbarm.set_label('I (K. km.$s^{-1}$)', size=10)
            mm = axm.imshow(maskm, cmap='gray', alpha=0.6)
        if i == 1:
            imm.set_clim(4., 9)
            cbarm.set_label('$V_{LSR}$ (K. km.$s^{-1}$)', size=10)
        if i == 2:
            imm.set_clim(0.1, 0.9)
            cbarm.set_label('$V_{LSR}$ ($km^2$.$s^{-2}$)', size=10)
        if i == 3:
            imm.set_clim(0.3, 2)
            cbarm.set_label('$\sigma _{V_{LSR}}$ ( km.$s^{-1}$)', size=10)

        plt.title(title[i])
        axm.invert_yaxis()  # Flips the Y axis

        axm.set_xlabel("RA (J2000)", fontsize=10)
        axm.set_ylabel("DEC (J2000)", fontsize=10)
        # Add a colorbar

        # Overplot column density contours
        cont_levels1m = [np.nanmax(momentm[0].hdu.data)*i for i in (0.3, 0.5, 0.7)]
        CS2m = axm.contour(momentm[0].hdu.data, levels=cont_levels1m, colors='black', linewidths=1)

        loopm += 1
        st = 1
        vl = moment_1m.hdu.data
        sg = linewidth_sigmam.hdu.data
        velm = []
        sigmam = []
    for y in range(0, ysm, st):
        for x in range(0, xsm, st):
            if maskm[y, x] == 3:
                sigmam.append(float(sg[y, x]))
                velm.append(float(vl[y, x]))
    print( 'mimir observations ---------------------------')
    print( "max/min of velocity", np.max(velm), np.min(velm))
    print( "average of velocity", np.average(velm))

    print( "std of velocity", np.nanstd(velm))
    print( "median of sigma", np.nanmedian(sigmam))
    print( "std of sigma", np.nanstd(sigmam))

    sig =[]
    vel = []
    for i in range(len(sigmam)):
        if sigmam[i]>0 and sigmam[i]<2:
            sig.append(sigmam[i])
            vel.append(velm[i])
    print( "std of velocity", np.nanstd(vel))
    print( "std of sigma", np.nanstd(sig))

    print( "IQR", iqr(sig))

    plt.figure()

    plt.hist(sig,bins=100)
    plt.xlabel("SIGMA (KM/S)")
    plt.ylabel("Pixels")


'''
need description
'''


def map_data(herschel_map, ra, dec, pa_s, sig_pa_s, deg_s, sig_deg_s, vlr, vlr_disp, Sig_linew):

    herschel = fits.open(get_pkg_data_filename(herschel_map))[0]
    f = fits.getdata(herschel_map, ext=0)
    f = np.log10(f)
    (ys, xs) = f.shape

    w = WCS(herschel.header)
    r, d = w.all_world2pix(np.array(ra), np.array(dec), 1)
    radec = []
    print( len(d))
    r = r.astype(int)
    d = d.astype(int)
    for i in range(len(d)):
        radec.append([int(r[i]), int(d[i])])
    #print radec
    print( "number of point ", len(radec))
    pa_so = np.zeros((ys, xs), dtype=np.float)
    sig_pa_so = np.zeros((ys, xs), dtype=np.float)
    deg_so = np.zeros((ys, xs), dtype=np.float)
    sig_deg_so = np.zeros((ys, xs), dtype=np.float)
    Vlr = np.zeros((ys, xs), dtype=np.float)
    vlr_d = np.zeros((ys, xs), dtype=np.float)
    sigma_linew = np.zeros((ys, xs), dtype=np.float)

    px_v, py_v = [], []
    for j in range(len(radec)):
        print( radec[j][0], radec[j][1], f[radec[j][1], radec[j][0]], pa_s[j], vlr[j], vlr_disp[j])
        pa_so[radec[j][1], radec[j][0]] = float(pa_s[j])
        sig_pa_so[radec[j][1], radec[j][0]] = float(sig_pa_s[j])

        deg_so[radec[j][1], radec[j][0]] = float(deg_s[j])
        sig_deg_so[radec[j][1], radec[j][0]] = float(sig_deg_s[j])

        Vlr[radec[j][1], radec[j][0]] = float(vlr[j])
        vlr_d[radec[j][1], radec[j][0]] = float(vlr_disp[j])
        sigma_linew[radec[j][1], radec[j][0]] = float(Sig_linew[j])
        if vlr[j] > 0:
            px_v.append(radec[j][0])
            py_v.append(radec[j][1])
    return px_v, py_v, Vlr, vlr_d, sigma_linew, pa_so, sig_pa_so, deg_so, sig_deg_so


px_v, py_v, Vlr, vlr_d, sigma_linew, pa_so, sig_pa_so, deg_so, sig_deg_so = map_data(
    herschel_map, ra, dec, pa_s, sig_pa_s, deg_s, sig_deg_s, vlr, vlr_disp, Sig_linew)
'''
Need description
these output will be used to derive the magnetic field in each region
px_v,py_v: are pixel position of the molecular line meseaurement that are also meseared in sofia data
Vlr: moment 1
vlr_d: moment_2
sigma_linew: sigma with line width
'''
if plot_all_maps == 1:
    '''1/. reprojection of planck and herschel '''
    herschel = fits.open(get_pkg_data_filename(herschel_map))[0]
    f = fits.getdata(herschel_map, ext=0)
    f = np.log10(f)
    print( '--------------- Taurus ------------------------------------------')
    pol_angle = planck_map
    # ------------------------------------------------------------------"
    print( 'processing to project the Planck map in Herschel')
    planck = fits.open(get_pkg_data_filename(pol_angle))[0]
    p_a, footprint = reproject_interp(planck, herschel.header)
    # -------------------------------------------------------------------
    print('Planck map are reporojected to Heschel map')

    #####################################
    p_a = p_a+np.pi/2
    ######################################
    #Operations in data
    print( 'Polarization data are rotated to get Magnetic field oriontation')
    (ys, xs) = p_a.shape
    print( "the vectors in map are magnetic field vectors obtained by rotating polarization data by 90 deg")
    print( "prepare sofia polarization data")
    st = 1
    r2 = 12
    fig = aplpy.FITSFigure(herschel, hdu=0, subplot=(1, 1, 1))
    linelist = []
    listhawc = []
    pa_sofia_all = []
    pa_planck_all = []
    f_all_r = []
    for y in range(0, ys, st):
        for x in range(0, xs, st):
            if pa_so[y, x] > 0 or pa_so[y, x] < 0:

                x1 = x+r2*np.cos(float(pa_so[y, x])*np.pi/180-np.pi/2)
                y1 = y+r2*np.sin(float(pa_so[y, x])*np.pi/180-np.pi/2)
                x2 = x-r2*np.cos(float(pa_so[y, x])*np.pi/180-np.pi/2)
                y2 = y-r2*np.sin(float(pa_so[y, x])*np.pi/180-np.pi/2)
                x_world, y_world = fig.pixel2world([x1, x2], [y1, y2])
                line = np.array([x_world, y_world])
                linelist.append(line)
                listhawc.append(line)
                f_all_r.append(c[y, x])
                pa_sofia_all.append(pa_so[y, x])
                pa_planck_all.append(p_a[y, x]*180/np.pi)
    print( "Sofia done !. preparing planck data ")
    st2 = 100
    r = 50
    pl = []
    linlist = []
    planck = []
    for y in range(0, ys, st2):
        for x in range(0, xs, st2):
            x1 = x+r*np.cos(float(p_a[y, x])-np.pi/2)
            y1 = y+r*np.sin(float(p_a[y, x])-np.pi/2)
            x2 = x-r*np.cos(float(p_a[y, x])-np.pi/2)
            y2 = y-r*np.sin(float(p_a[y, x])-np.pi/2)
            x_world, y_world = fig.pixel2world([x1, x2], [y1, y2])
            lin = np.array([x_world, y_world])
            linlist.append(lin)
            planck.append(p_a[y, x]*180/np.pi)

    #print "std planck sub-region ", np.std(planck_pa_sub),planck_pa_sub
    #print 'max of density in all sofia measurements',np.nanmax(f_all_r)
    a = fig.show_contour(herschel, levels=[3*10**21, 6.7*10 **
                                           21, 10**22, 10**23], overlap=True, alpha=0.5)
    s = fig.show_lines(linelist, layer='vectdors', color='black', alpha=0.8, linewidth=1)
    pl = fig.show_lines(linlist, layer='vect', color='yellow', alpha=0.5, linewidth=1.3)

    sov = linelist
    im = plt.imshow(f, cmap='gist_stern')
    d = fig.recenter(64.6, 27.45, width=0.43, height=0.43)
    im.set_clim(20.5, 22.5)
    plt.scatter(px_v, py_v, s=20, c='red', marker='o', alpha=0.6)
    fig.add_scalebar(0.041, color='white', linewidth=1.5,
                     weight='bold', corner='bottom right')  # length has to be specified
    fig.scalebar.set_label('0.1 pc')
    # Overplot column density contours

    cax = plt.axes([0.86, 0.11, 0.05, 0.77])
    fig.add_beam(label='SOFIA beam')
    fig.beam.set_color('white')
    fig.beam.set_hatch('+')
    plt.colorbar(mappable=im,cax = cax)

plt.show()
def DCF(vel_dis, pa_d, n_H2):
    C = np.sqrt(8*np.log(2))
    #print 'constant',C
    B = 9.3*np.sqrt(n_H2)*(vel_dis/pa_d)*C
    return B


if plot_velocity == 1:
    'Read file contaning the result using gildas '
    VLR1 = np.loadtxt('result_fit_REG1_all.txt', comments='!')[:, 6]
    VLR2 = np.loadtxt('result_fit_REG2.txt', comments='!')[:, 6]
    VLR3 = np.loadtxt('result_fit_REG3.txt', comments='!')[:, 6]
    VLR4 = np.loadtxt('result_fit_r4.txt', comments='!')[:, 6]

    dis_VLR1 = np.array((np.loadtxt('result_fit_REG1_all.txt',
                                    comments='!')[:, 8]))/(np.sqrt(8*log(2)))
    dis_VLR2 = np.array((np.loadtxt('result_fit_REG2.txt', comments='!')[:, 8]))/(np.sqrt(8*log(2)))
    dis_VLR3 = np.array(np.loadtxt('result_fit_r3.txt', comments='!')[:, 8])/(np.sqrt(8*log(2)))
    dis_VLR4 = np.array(np.loadtxt('result_fit_r3.txt', comments='!')[:, 8])/(np.sqrt(8*log(2)))

    'This is to plot magnetic field angles for Planck and Hawc+ '
    fig3 = plt.figure(figsize=(12, 8))
    to_plot = [VLR1, VLR2, VLR3, VLR4, dis_VLR1, dis_VLR2, dis_VLR3, dis_VLR4]

    print( 'number of mesureaments')
    print( 'region1', len(VLR1), len(dis_VLR1))
    print( 'region2', len(VLR2), len(dis_VLR2))
    print( 'region3', len(VLR3), len(dis_VLR3))
    print( 'region4', len(VLR4), len(dis_VLR4))
    title = ['Region 1', 'Region 2', 'Region 3', 'Region 4']
    loop = 1
    num_bins = 15
    row, col = 2, 4
    l = 0
    for j in range(8):

        axes = fig3.add_subplot(row, col, loop)
        # the histogram of the data
        axes.hist(to_plot[j], num_bins, alpha=0.5)
        axes.set_xlabel('Velocity in km/s')
        binwidth = 0.1
        #num_bins = np.int((np.max(to_plot[j])-np.min(to_plot[j]))/binwidth)
        num_bins = 8
        axes.hist(to_plot[j], num_bins, alpha=0.5, color='gray', edgecolor='red', hatch='x2')
        if l < 4:
            plt.title(title[l])
        if j > 3:

            axes.set_xlabel('Velocity dispersion (km/s)')

        if j == 0 or j == 4:
            axes.set_ylabel('number of gaussian fits')
        l = l+1
        loop += 1
    v_co1 = []
    v_co2 = []
    v_co3 = []
    d_co1 = []
    d_co2 = []
    d_co3 = []
    l = 0
    st = 1

    for l in range(0, len(VLR3)-3, st):

        v_co1.append(VLR3[l])
        v_co2.append(VLR3[l+1])
        v_co3.append(VLR3[l+2])
        d_co1.append(dis_VLR3[l])
        d_co2.append(dis_VLR3[l+1])
        d_co3.append(dis_VLR3[l+2])
        print( VLR2[l], dis_VLR3[l])
        #print VLR2[l+1]
        #print VLR2[l+2]
    vel_c1 = []
    vel_c2 = []
    vel_c1d = []
    vel_c2d = []
    vel_c3 = []
    vel_c3d = []
    for k in range(0, len(VLR3), st):
        if int(VLR3[k]) == 4:
            vel_c1.append(VLR3[k])
            vel_c1d.append(dis_VLR3[k])
            print( dis_VLR3[k])
        if int(VLR2[k]) == 5:
            vel_c2.append(VLR3[k])
            vel_c2d.append(dis_VLR3[k])
        if int(VLR2[k]) == 6:
            vel_c3.append(VLR3[k])
            vel_c3d.append(dis_VLR3[k])

    print( "first component mean+/-std", np.average(v_co1), np.std(v_co1), np.average(d_co1))
    print( "second component mean+/-std", np.average(v_co2), np.std(v_co2), np.average(d_co2))
    print( "third component mean+/-std", np.average(v_co3), np.std(v_co3), np.average(d_co3))
    print( "----------------- First component", np.average(
        vel_c1), np.std(vel_c1), np.average(vel_c1d))
    print("-----------------Second component", np.average(
        vel_c2), np.std(vel_c2), np.average(vel_c2d))
    print( "-----------------3 component", np.average(vel_c3), np.std(vel_c3), np.average(vel_c3d))
    print( 'number of velocity mesurement ', len(VLR3), len(vel_c1), len(vel_c2), len(vel_c3))


def DFC_sub_region(region, Vlr, vlr_d, sigma_linew, pa_so, sig_pa_so, deg_so, sig_deg_so, VLR, dis_VLR):

    #mask6 = np.zeros((ys, xs), dtype=np.float)
    w = WCS(herschel.header)
    radec = []
    for i in range(4):
        r, d = w.all_world2pix(region[i][0], region[i][1], 1)
        radec.append([int(r), int(d)])

    print( "pixel positions of the region in herschel map", radec)
    mask = np.zeros((ys, xs), dtype=np.float)
    polys = np.array([radec])
    import cv2
    import os

    for i in polys.tolist():
        cv2.fillConvexPoly(mask, np.array(i), 1)
    mask[mask == 0] = np.nan
    st = 1
    r2 = 4
    linelist = []
    den = []
    pa_disp = []
    V = []
    Sig = []
    Vd = []
    planck_Bf = []
    hawc = []
    f_d = []
    # disV: is the disperssion of velocity obtained with 30m data
    print( "density,    vlr(km/s)   , Velocity_disperssion(km/s)  , p_a(deg), p_d(per)"
)
    for y in range(0, ys, st):
        for x in range(0, xs, st):
            # if (mask3[y,x]==1 or mask3[y,x]==1 or mask4[y,x]==1 or mask6[y,x]==1) and (pa_so[y,x]>0 or pa_so[y,x]<0) :
            # if (mask[y,x]==1) and Vlr[y,x]>0 and vlr_d[y,x]>0 and pa_so[y,x]/sig_pa_so[y,x]<3 and f[y,x]<density_cr:
            # if (mask[y,x]==1) and Vlr[y,x]>0 and vlr_d[y,x]>0 and f[y,x]<density_cr and pa_so[y,x]/sig_pa_so[y,x]>=3 :
            if Vlr[y, x] > 0 and vlr_d[y, x] > 0:
                pa_disp.append(pa_so[y, x])
                V.append(Vlr[y, x])
                Vd.append(vlr_d[y, x])
                Sig.append(sigma_linew[y, x])
                den.append(f[y, x])
                planck_Bf.append(p_a[y, x]*180/np.pi)
                hawc.append(pa_so[y, x]-90)
                f_d.append(f[y, x])
                print( f[y, x], Vlr[y, x], vlr_d[y, x], sigma_linew[y, x], pa_so[y, x], deg_so[y, x])

    print( 'mean column density', np.nanmean(den), "std", np.std(
        den), "number of meseaurement in this region ", len(den))
    print( 'polarization disperssion', format(np.std(pa_disp), '.2g'), format(np.std(hawc), '.2g'))
    print( pa_disp)
    print('mean of centeroid velocity', format(np.average(V), '.2g'))
    print( 'mean of density', format(np.average(f_d), '.2g'))
    print( 'velocity disperssion / average of sqrt(moment2)', format(np.nanmean(Vd), '.2g'), "standard deviation of velocity ", format(
        np.std(V), '.2g'), "average of sima using line width ", format(np.nanmean(Sig), '.2g'))
    print( "NUMBER OF MEASUREMENTS #######", "IRAM30 gaussian fits ", len(
        VLR), len(dis_VLR), "SOFIA HAWC", len(den))
    print('number of IRAM measurements ', len(V), len(Vd))
    print( "Methode 1 : using std of the velocity value in this region")

    vel_dis = np.std(V)
    pa_d = np.std(pa_disp)
    B = DCF(vel_dis, pa_d, n_H2)

    print( "Magnetic field strength in this region is : ", format(B, '.2g'), "uG")

    print( "Methode 2 : using variance of the velocity (M2)")
    vel_dis = np.average(Vd)
    pa_d = np.std(pa_disp)
    B = DCF(vel_dis, pa_d, n_H2)

    print( "Magnetic field strength in this region is : ", format(B, '.2g'), "uG")

    print( "Methode 3 : using using line width")
    vel_dis = np.average(Sig)
    pa_d = np.std(pa_disp)
    B = DCF(vel_dis, pa_d, n_H2)

    print( "Magnetic field strength in this region is : ", format(B, '.2g'), "uG")

    print( "Methode 4 : using Hacar et al. 2013 data ")
    vel_dis = 0.2  # Hacar data
    pa_d = np.std(pa_disp)
    B = DCF(vel_dis, pa_d, n_H2)
    print( "Magnetic field strength in this region is : ", format(B, '.2g'), "uG")

    print( "Methode 5 : GILDAS ANALYSIS , np.std(VLR) ")
    vel_dis = np.std(VLR)  # Hacar data
    pa_d = np.std(pa_disp)
    B = DCF(vel_dis, pa_d, n_H2)
    print( "Magnetic field strength in this region is : ", format(B, '.2g'), "uG")
    print( "Methode 6 : GILDAS ANALYSIS , np.average(disp_VLR) ")
    vel_dis = np.average(dis_VLR3)  # Hacar data
    pa_d = np.std(pa_disp)
    B = DCF(vel_dis, pa_d, n_H2)
    print( "Magnetic field strength in this region is : ", format(B, '.2g'), "uG")
    return B, V, planck_Bf, hawc


B, V, planck_Bf, hawc = DFC_sub_region(
    region, Vlr, vlr_d, sigma_linew, pa_so, sig_pa_so, deg_so, sig_deg_so, VLR3, dis_VLR3)
if plot_planck == 1:
    fig = aplpy.FITSFigure(herschel, hdu=0, subplot=(1, 1, 1))
    density = []
    sk_number = []
    linelist = []
    linlist = []
    mask = np.zeros((ys, xs), dtype=np.int8)
    a, b = 0, 0
    polys = np.array([[[2340, 2722], [4435, 4078], [5961, 3585], [3280, 1690]]])
    import cv2
    for i in polys.tolist():
        cv2.fillConvexPoly(mask, np.array(i), 1)
    a = fig.recenter(64.8278728, 27.3861307, width=3, height=3)
    # reg=[[3924,3794],[5127,3794],[5127,2597],[3924,2605]]
    # reg=[[2860,3558],[3962,4510],[5599,3419],[3930,1953]]
    reg = [[2720, 3258], [4561, 5056], [6305, 2980], [4390, 1172]]
    radec = [[3159, 2862], [4614, 4189], [5866, 3076], [4261, 1536]]
    (ys, xs) = herschel.data.shape
    mask = np.zeros((ys, xs), dtype=np.float)
    polys = np.array([reg])
    for i in range(4):
        import cv2
        import os

        for i in polys.tolist():
            cv2.fillConvexPoly(mask, np.array(i), 1)
        mask[mask == 0] = np.nan
    st = 200
    r2 = 150
    planck = []
    #####################
    regp = [[3197, 2898], [4768, 4192], [5181, 4157], [5730, 3268],
            [3936, 1881]]  # premiere estimation du champ magnetic
    #regp= [[2512,3548],[4198,4849],[5707,3199],[3916,1892]]
    #regp= [[2512,3548],[4148,4859],[5663,3241],[3878,1957]]

    (ys, xs) = herschel.data.shape
    maskp = np.zeros((ys, xs), dtype=np.float)
    polysp = np.array([regp])
    for i in range(4):
        import cv2
        import os
        for i in polysp.tolist():
            cv2.fillConvexPoly(maskp, np.array(i), 1)
        maskp[maskp == 0] = np.nan

    ##############################
    st = 1000
    dens = []

    w = WCS(herschel.header)
    radecm = []
    for i in range(4):
        rm, dm = w.all_world2pix(regmimir[i][0], regmimir[i][1], 1)
        radecm.append([int(rm), int(dm)])
        maskm = np.zeros((ys, xs), dtype=np.float)
        polysm = np.array([radecm])
        import cv2
        import os

    for im in polysm.tolist():
        cv2.fillConvexPoly(maskm, np.array(im), 3)
        maskm[maskm == 0] = np.nan
    st = 200
    r2 = 100
    planckp = []
    planck = []
    pl = []
    density_p = []
    for y in range(0, ys, st):
        for x in range(0, xs, st):
            if maskp[y, x] == 1:
                x1 = x+r2*np.cos(float(p_a[y, x])-np.pi/2)
                y1 = y+r2*np.sin(float(p_a[y, x])-np.pi/2)
                x2 = x-r2*np.cos(float(p_a[y, x])-np.pi/2)
                y2 = y-r2*np.sin(float(p_a[y, x])-np.pi/2)
                x_world, y_world = fig.pixel2world([x1, x2], [y1, y2])
                line = np.array([x_world, y_world])
                linelist.append(line)
                planckp.append(p_a[y, x]*180/np.pi)

                #pl.append(float(p_a[y, x]))
    linlist = []
    lin = []

    st = 200
    r2 = 100
    x, y = 0, 0

    for y in range(0, ys, st):
        for x in range(0, xs, st):
            if maskm[y, x] == 3:
                x1 = x+r2*np.cos(float(p_a[y, x])-np.pi/2)
                y1 = y+r2*np.sin(float(p_a[y, x])-np.pi/2)
                x2 = x-r2*np.cos(float(p_a[y, x])-np.pi/2)
                y2 = y-r2*np.sin(float(p_a[y, x])-np.pi/2)
                x_world, y_world = fig.pixel2world([x1, x2], [y1, y2])
                lin = np.array([x_world, y_world])
                linlist.append(lin)
                #planck.append(p_a[y, x]*180/np.pi)
                pl.append(float(p_a[y, x]))
                density_p.append(c[y, x])
                planck.append(p_a[y, x]*180/np.pi)

    #sap = fig.show_lines(listhawc, layer='vectorshawc', color='orange', alpha=1, linewidth=1)
    s = fig.show_lines(linelist, layer='vectorsp', color='orange', alpha=1, linewidth=1)
    #a=fig.show_contour(herschel,levels=[(21)*i for i in (0.5,0.7)],overlap=True,color='black')
    a = fig.show_contour(herschel, levels=[3*10**21, 6.7*10 **
                                           21, 10**22, 10**23], overlap=True, alpha=0.5)
    sigf = fig.show_lines(linlist, layer='vectorsplanckestim',
                          color='yellow', alpha=1, linewidth=1.2)

    im = plt.imshow(f, cmap='gist_stern')
    maskp[maskp == 1] = 23
    # ii=plt.imshow(maskm,cmap='gray',alpha=0.6)
    cax = plt.axes([0.86, 0.11, 0.05, 0.77])
    plt.colorbar(mappable=im, cax=cax)
    im.set_clim(20.5, 22.5)
    print( "std Planck region --------------", np.average(planck), np.nanstd(planck))
    #print "std Planck big region ", np.nanstd(pl), len(pl)
    print( "average density in yellow vectors", np.average(
        np.array(density_p)/1.e+21)*1.e+21, np.std(np.array(density_p)/.1e+21), len(planck))

    plt.scatter(px_v, py_v, s=40, c='red', marker='o', alpha=0.6)
    fig.add_scalebar(0.164, color='black', linewidth=1.5,
                     weight='bold', corner='bottom right')  # length has to be specified
    fig.scalebar.set_label('0.4 pc')
    # Overplot column density contours

    cax = plt.axes([0.86, 0.11, 0.05, 0.77])
    fig.add_beam(label='SOFIA beam')
    fig.beam.set_color('black')
    fig.beam.set_hatch('+')
if plot_stat_B == 1:
    'This is to plot magnetic field angles for Planck and Hawc+ '
    fig3 = plt.figure(figsize=(12, 8))
    print( 'number of point in mask', len(hawc), len(planck_Bf))
    print( 'all', len(pa_planck_all), len(pa_sofia_all))
    title = ['B-field from HAWC+ data in this region ', 'B-field from planck in this region',
             'velocity in this region ', "B-field angle with Sofia in all regions", "B-field angle with Planck in all regions"]
    loop = 1
    #num_bins = 5
    binwidth = 15.
    row, col = 1, 3
    num_bins_sofia = np.int((np.max(pa_sofia_all)-np.min(pa_sofia_all))/binwidth)
    axes = fig3.add_subplot(row, col, 1)
    # the histogram of the data
    entries, edges, _ = plt.hist(pa_sofia_all, bins=num_bins_sofia,
                                 alpha=0.7, color='yellow', edgecolor='red', hatch='x2')
    # calculate bin centers

    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    plt.errorbar(bin_centers, entries, yerr=np.sqrt(entries), fmt='b.')
    print( 'number of sofia mesearements ', len(pa_sofia_all))
    print( "bin width histogramme sofia", num_bins_sofia)
    axes.set_xlabel('B-field angles')
    axes.set_ylabel('Number of independent measurements')
    plt.xlim(-90, 90)
    #axes.set_title("B-field angle obtained with HAWC+ observations")

    axes = fig3.add_subplot(row, col, 2)
    pla = []
    for i in range(len(pl)):
        if pl[i] > np.pi:
            pla.append(np.pi - pl[i])
        else:
            pla.append(pl[i])
    pl = pla
    print( len(pla)*180/np.pi, np.mean(pla)*180/np.pi, np.std(pla)*180/np.pi, np.max(pla)*180/np.pi)
    num_bins_planck = np.int((np.max(np.array(pl)*180/np.pi) -
                              np.min(np.array(pl)*180/np.pi))/binwidth)
    plt.axvline(x=-81, dash_joinstyle='round', c='red', linestyle='-.', linewidth=2)
    plt.text(-74, 45, 'Filament', rotation=90, color='red', weight='bold')
    entries, edges, _ = plt.hist(np.array(pl)*180/np.pi, bins=num_bins_planck,
                                 alpha=0.7, color='green', edgecolor='red', hatch='x2')
    # calculate bin centers
    print( "Number of value per bin", entries)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    plt.errorbar(bin_centers, entries, yerr=np.sqrt(entries), fmt='g.')
    axes.set_xlabel('B-field angles')
    plt.xlim(-90, 180)
    axes.set_ylabel('Number of independent measurements')
    plt.axvline(x=-81, dash_joinstyle='round', c='red', linestyle='-.', linewidth=2)
    plt.text(-74, 37, 'Filament', rotation=90, color='red', weight='bold')
    print( 'number of Planck mesearements ', len(pl), np.min(pl), np.max(pl))
    print( "bin width histogramme Planck", num_bins_planck)
    axes = fig3.add_subplot(row, col, 3)
    #axes.hist(np.array(planck)*180/np.pi,num_bins_planck,alpha=0.9, color = 'green',edgecolor = 'red', hatch = 'x2',label='Planck')
    ntries, edges, _ = plt.hist(np.array(pl)*180/np.pi, num_bins_planck,
                                alpha=1, color='green', edgecolor='red', hatch='x2', label='Planck')
    # calculate bin centers

    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    plt.errorbar(bin_centers, entries, yerr=np.sqrt(entries), fmt='g.')
    entries, edges, _ = plt.hist(pa_sofia_all, bins=num_bins_sofia, alpha=0.5,
                                 color='yellow', edgecolor='red', hatch='x2', label='HAWC+')
    # calculate bin centers
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    plt.errorbar(bin_centers, entries, yerr=np.sqrt(entries), fmt='b.')
    #axes.hist(pa_sofia_all,num_bins_sofia,alpha=0.5, color = 'yellow',edgecolor = 'red', hatch = 'x2',label='HAWC+')
    plt.axvline(x=-81, dash_joinstyle='round', c='red', linestyle='-.', linewidth=2)
    plt.text(-74, 45, 'Filament', rotation=90, color='red', weight='bold')
    plt.legend()
    axes.set_xlabel('B-field angles')
    axes.set_ylabel('Number of independent measurements')
    plt.xlim(-90, 90)

    #axes.set_title("B-field angle obtained with Planck  and HAWC+ observations")


plt.show()
