if plot_stat_B == 1:
    'This is to plot magnetic field angles for Planck and Hawc+ '
    fig3 = plt.figure(figsize=(12, 8))
    print 'number of point in mask', len(hawc), len(planck_Bf)
    print 'all', len(pa_planck_all), len(pa_sofia_all)
    title = ['B-field from HAWC+ data in this region ', 'B-field from planck in this region',
             'velocity in this region ', "B-field angle with Sofia in all regions", "B-field angle with Planck in all regions"]
    loop = 1
    #num_bins = 5
    binwidth = 15.
    row, col = 1, 3
    num_bins_sofia = np.int((np.max(pa_sofia_all)-np.min(pa_sofia_all))/binwidth)
    axes = fig3.add_subplot(row, col, 1)
    # the histogram of the data
    entries, edges, _ = plt.hist(pa_sofia_all, bins=num_bins_sofia,
                                 alpha=0.7, color='yellow', edgecolor='red', hatch='x2')
    # calculate bin centers

    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    plt.errorbar(bin_centers, entries, yerr=np.sqrt(entries), fmt='b.')
    print 'number of sofia mesearements ', len(pa_sofia_all)
    print "bin width histogramme sofia", num_bins_sofia
    axes.set_xlabel('B-field angles')
    axes.set_ylabel('Number of independent measurements')
    plt.xlim(-90, 90)
    #axes.set_title("B-field angle obtained with HAWC+ observations")

    axes = fig3.add_subplot(row, col, 2)
    pla = []
    for i in range(len(pl)):
        if pl[i] > np.pi:
            pla.append(np.pi - pl[i])
        else:
            pla.append(pl[i])
    pl = pla
    print len(pla)*180/np.pi, np.mean(pla)*180/np.pi, np.std(pla)*180/np.pi, np.max(pla)*180/np.pi
    num_bins_planck = np.int((np.max(np.array(pl)*180/np.pi) -
                              np.min(np.array(pl)*180/np.pi))/binwidth)
    plt.axvline(x=-81, dash_joinstyle='round', c='red', linestyle='-.', linewidth=2)
    plt.text(-74, 45, 'Filament', rotation=90, color='red', weight='bold')
    entries, edges, _ = plt.hist(np.array(pl)*180/np.pi, bins=num_bins_planck,
                                 alpha=0.7, color='green', edgecolor='red', hatch='x2')
    # calculate bin centers
    print "Number of value per bin", entries
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    plt.errorbar(bin_centers, entries, yerr=np.sqrt(entries), fmt='g.')
    axes.set_xlabel('B-field angles')
    plt.xlim(-90, 180)
    axes.set_ylabel('Number of independent measurements')
    plt.axvline(x=-81, dash_joinstyle='round', c='red', linestyle='-.', linewidth=2)
    plt.text(-74, 37, 'Filament', rotation=90, color='red', weight='bold')
    print 'number of Planck mesearements ', len(pl), np.min(pl), np.max(pl)
    print "bin width histogramme Planck", num_bins_planck
    axes = fig3.add_subplot(row, col, 3)
    #axes.hist(np.array(planck)*180/np.pi,num_bins_planck,alpha=0.9, color = 'green',edgecolor = 'red', hatch = 'x2',label='Planck')
    ntries, edges, _ = plt.hist(np.array(pl)*180/np.pi, num_bins_planck,
                                alpha=1, color='green', edgecolor='red', hatch='x2', label='Planck')
    # calculate bin centers

    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    plt.errorbar(bin_centers, entries, yerr=np.sqrt(entries), fmt='g.')
    entries, edges, _ = plt.hist(pa_sofia_all, bins=num_bins_sofia, alpha=0.5,
                                 color='yellow', edgecolor='red', hatch='x2', label='HAWC+')
    # calculate bin centers
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    plt.errorbar(bin_centers, entries, yerr=np.sqrt(entries), fmt='b.')
    #axes.hist(pa_sofia_all,num_bins_sofia,alpha=0.5, color = 'yellow',edgecolor = 'red', hatch = 'x2',label='HAWC+')
    plt.axvline(x=-81, dash_joinstyle='round', c='red', linestyle='-.', linewidth=2)
    plt.text(-74, 45, 'Filament', rotation=90, color='red', weight='bold')
    plt.legend()
    axes.set_xlabel('B-field angles')
    axes.set_ylabel('Number of independent measurements')
    plt.xlim(-90, 90)
