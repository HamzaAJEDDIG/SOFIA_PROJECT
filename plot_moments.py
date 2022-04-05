from Used_modules import *
from manip_SOFIA_DATA import *
from Input import *
from manip_CUBE_DATA import *

cube, sub_cub, moment_0, moment_1, moment_2, linewidth_sigma = extract_moments(cub)

x, y, ra, dec, pa_s, sig_pa_s, deg_s, sig_deg_s, vlr, vlr_disp, Sig_linew = sofia_dat(
    table, moment_1, moment_2, linewidth_sigma)

def plot_moments():
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