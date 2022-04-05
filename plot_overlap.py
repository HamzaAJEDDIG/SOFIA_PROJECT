from Used_modules import *
from manip_SOFIA_DATA import *
from Input import *
from manip_CUBE_DATA import *

cube, sub_cub, moment_0, moment_1, moment_2, linewidth_sigma = extract_moments(cub)

x, y, ra, dec, pa_s, sig_pa_s, deg_s, sig_deg_s, vlr, vlr_disp, Sig_linew = sofia_dat(
    table, moment_1, moment_2, linewidth_sigma)

def plot_moment0():
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
    #plt.show()