from Used_modules import *
from Input import *
def plot_moment_mimir():
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
