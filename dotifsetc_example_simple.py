#python3 dotifsetc_example_simple.py
from dotifsetc import dotifsetc

#run DOTIFS ETC with below conditions:
#Exposure time: 900 sec
#Source: Sc type galaxy, r band surface brightness 17 mag/arcsec2
#Source redshift: z=0.04
#Sky: Halfmoon, r band surface brightness 22 mag/arcsec2
#Wavelength bin size: 2.5 pixels (=3.08 Angstrom)
res=dotifsetc(exptime=900, source='obj_sc', z=0.04,
		magnitude=17, skymagnitude=22, stype=1,
		pixel=2.5,
		oname='dotifs_etc.jpg', save=True)

#change exposure time to 3600 sec
#change output name as '3600.ps'
#change band to 'g' band (both source and sky)
#change plot range: from 5000 to 7000 Angstrom
#turn off sky subtraction
#change sky g band magnitude to 18 
res.exptime=3600
res.oname='3600.jpg'
res.band='g'
res.skymagnitude=18
res.skysub=False
res.run()
