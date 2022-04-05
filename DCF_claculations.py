from Input import *
def DCF(vel_dis, pa_d, n_H2):
    C = np.sqrt(8*np.log(2))
    #print 'constant',C
    B = 9.3*np.sqrt(n_H2)*(vel_dis/pa_d)*C
    return B


def plot_velocity():
    'Read file contaning the result using gildas '
    VLR1 = np.loadtxt('./data/result_fit_REG1_all.txt', comments='!')[:, 6]
    VLR2 = np.loadtxt('./data/result_fit_REG2.txt', comments='!')[:, 6]
    VLR3 = np.loadtxt('./data/result_fit_REG3.txt', comments='!')[:, 6]
    VLR4 = np.loadtxt('./data/result_fit_r4.txt', comments='!')[:, 6]

    dis_VLR1 = np.array((np.loadtxt('./data/result_fit_REG1_all.txt',
                                    comments='!')[:, 8]))/(np.sqrt(8*log(2)))
    dis_VLR2 = np.array((np.loadtxt('./data/result_fit_REG2.txt', comments='!')[:, 8]))/(np.sqrt(8*log(2)))
    dis_VLR3 = np.array(np.loadtxt('./data/result_fit_r3.txt', comments='!')[:, 8])/(np.sqrt(8*log(2)))
    dis_VLR4 = np.array(np.loadtxt('./data/result_fit_r3.txt', comments='!')[:, 8])/(np.sqrt(8*log(2)))

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

    return  VLR3, dis_VLR3
def DFC_sub_region(region, Vlr, vlr_d, sigma_linew, pa_so, sig_pa_so, deg_so, sig_deg_so, VLR, dis_VLR,p_a,VLR3, dis_VLR3):

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

    #print( 'mean column density', np.nanmean(den), "std", np.std(
    #    den), "number of meseaurement in this region ", len(den))
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

