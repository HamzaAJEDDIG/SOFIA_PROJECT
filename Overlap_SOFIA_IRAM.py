def sofia_dat(table,moment_1,moment_2,linewidth_sigma):
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
