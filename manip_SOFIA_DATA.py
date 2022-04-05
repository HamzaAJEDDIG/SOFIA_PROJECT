from Used_modules import *
from Input import *
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
