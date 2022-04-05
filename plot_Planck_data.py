from Used_modules import *
from Input import *
def plot_planck(p_a,px_v, py_v):
    fig = aplpy.FITSFigure(herschel, hdu=0, subplot=(1, 1, 1))
    density = []
    sk_number = []
    linelist = []
    linlist = []
    (ys, xs) = f.shape
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