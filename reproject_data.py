from Input import *
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

