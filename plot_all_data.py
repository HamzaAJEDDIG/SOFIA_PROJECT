from Used_modules import *
from Input import *
def plot_all_maps(px_v, py_v, Vlr, vlr_d, sigma_linew, pa_so, sig_pa_so, deg_so, sig_deg_so):
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

    return px_v, py_v, Vlr, vlr_d, sigma_linew, pa_so, sig_pa_so, deg_so, sig_deg_so,p_a