# SOFIA_PROJECT
# Envirement and aims 
  - Code was built in python 3.8, all needed modules are in venv directory
  - This project combine different data sets from different telescopes such as Herschel, IRAM 30m, MIMIR instrument and Planck
  - The pipeline allows to compute the magnetic field stregth using DCF method                   #
                               the input :                                                                   
                                          1/. fits file of cub data containing velocity                      
                                          2/. fits file conataining the backgroud map                        
                                            (in this case Herscel data)                                      
                                          3/. table containing polarization data organize                    
                                                RA, DEC ,PA,  SIF_PA,  PDEG,  SIG_PDEG                       
                                          4/. fits file containing Planck data to be displayed (PA)          

# Results of this project are published in MINRAS 
# Authors of this code
# Hamza AJEDDIG 
# Authors of the paper 
# Li, Pak Shing  search by orcid ;  Lopez-Rodriguez, Enrique ;  Ajeddig, Hamza ;  André, Philippe ; McKee, Christopher F.  search by orcid ;  Rho, Jeonghee ;  Klein, Richard I.

# Abstract
Optical and infrared polarization mapping and recent Planck observations of the filametary cloud L1495 in Taurus show that the large-scale magnetic field is approximately perpendicular to the long axis of the cloud. We use the HAWC + polarimeter on SOFIA to probe the complex magnetic field in the B211 part of the cloud. Our results reveal a dispersion of polarization angles of 36°, about five times that measured on a larger scale by Planck. Applying the Davis-Chandrasekhar-Fermi (DCF) method with velocity information obtained from Institut de Radioastronomie Millimétrique 30 m C18O(1-0) observations, we find two distinct sub-regions with magnetic field strengths differing by more than a factor 3. The quieter sub-region is magnetically critical and sub-Alfvenic; the field is comparable to the average field measured in molecular clumps based on Zeeman observations. The more chaotic, super-Alfvenic sub-region shows at least three velocity components, indicating interaction among multiple substructures. Its field is much less than the average Zeeman field in molecular clumps, suggesting that the DCF value of the field there may be an underestimate. Numerical simulation of filamentary cloud formation shows that filamentary substructures can strongly perturb the magnetic field. DCF and true field values in the simulation are compared. Pre-stellar cores are observed in B211 and are seen in our simulation. The appendices give a derivation of the standard DCF method that allows for a dispersion in polarization angles that is not small, present an alternate derivation of the structure function version of the DCF method, and treat fragmentation of filaments.




