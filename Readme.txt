################################################################################
The project aim is to provide data structure and analysis of SOFIA project about
magnetic field in star-forming regions B211.
It allows the user to project different data sets from various telescopes :
   - The IRAM 30m data @
   - SOFIA polarization data
   - Planck polarization data
   - Herchel column density Data
   - NIKA2-Pol polarization data (optional)

the input :
  1/. fits file of cub data containing velocity
  2/. fits file conataining the backgroud map
   (in this case Herscel data)
  3/. table containing polarization data organize
    (Coordinates RA-DEC, Eolarization angle, Error in the polarization angles, Polarization degree, Error on the polarization degree)
    RA, DEC ,PA,  SIF_PA,  PDEG,  SIG_PDEG
 4/. fits file containing Planck data to be displayed (PA)

################################################################################
