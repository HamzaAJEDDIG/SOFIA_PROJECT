from Used_modules import *
from manip_CUBE_DATA import *
from Input import *
def plot_spectra(cube, sub_cub, moment_0, moment_1, moment_2, linewidth_sigma,x, y, ra, dec, pa_s, sig_pa_s, deg_s, sig_deg_s, vlr, vlr_disp, Sig_linew):
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