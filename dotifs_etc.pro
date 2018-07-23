;+
;NAME:
;	dotifs_etc, exptime=exptime, band=band, magnitude=magnitude, $
;	oname=oname, galtemp=galtemp, z=z, inputflux=inputflux, inputwave=inputwave, bbtemp=bbtemp, calibration=calibration, $
;	stype=stype, skymagnitude=skymagnitude, $
;	ftype=ftype, asahi=asahi, $
;	noplot=noplot, plotrange=plotrange, ptype=ptype, wstep=wstep, pixel=pixel, cspline=cspline, wavearr=wavearr, scflag=scflag$
;	waveout=waveout, snrout=snrout, signalarr=signalarr, skysignalarr=skysignalarr, noisearr=noisearr, noiseskyarr=noiseskyarr, $
;	dtypes=dtypes, ofile=ofile, $
;	help=help
;
;-
;
; PURPOSE:
;       Calculate expected number of spectral count from model source with predicted DOTIFS throughput.
;	Details can be found in 'DOTIFS_SNR_calculator_document.pdf'
;
; INPUTS: 
;	exptime - exposure time in seconds
;		default: 900 seconds
;	band - photometric band which will be used to calculate source flux. Only SDSS ugriz bands are supported.
;		default: 'r'
;	magnitude - source magnitude at the defined photometric band. The template flux will be scaled based on this
;		input magnitude. (AB magnitude)
;		default: 17 mag
;
; OPTIONAL INPUTS: 
;	oname - name of output ps file
;		default: 'snr.ps'
;
;	galtemp - source galaxy spectral model template type. Available templates can be found under ./kc96 directory.
;		ETC looking for template file with a name of galtemp+'_template.ascii'
;		default: 'sc' 
;	z - redshift of the galtemp source. It is only affected when a source is selected by galtemp
;		default: 0.000
;
;	inputflux - User can use user defined source flux as a source of ETC in unit of erg/s/cm^2/Ang. 
;		inputwave should be provided as well.
;	
;	inputwave - wavelength array for inputflux in unit of angstrom
;
;	bbtemp - blackbody temperature in kelvin. If this input is given then ETC assume blackbody object at given temperature
;		as a source. galtemp or inputflux will be ignored.
;
;	stype - choose sky type index. Currently, Four sky template spectra are supported as below
;		['sky150701_newmoon_alt90.dat','sky160311_newmoon_alt45.dat','sky160311_halfmoon.dat','sky160311_fullmoon.dat']
;		Their index is from 0 to 3.
;		default: 0 (newmoon)
;
;	skymagnitude - if this input is given then model sky flux will be scaled on this magnitude at selected photometric band. Otherwise,
;		absolute flux defined by stype will be used.
;
;	plotrange - determine wavelength range of the output. This should be provided in array with two elements.
;		default: [3700.d, 7400.d]
;
;	ptype - User can use pre-defined plot range by ptype. There are three plot ranges as below.
;		1: [3700,7400], 2: [3700,3760], 3: [7100,7400]
;		if thie input is given then plotrange will be ignored.
;	
;	ftype - Choose filter type. ETC assumes three diffrent ideal filter cut-on range at the blue side as below. filter transmission
;		will go up from 0.001 to 0.9 in this wavelength range. It affects second order contamination count.
;		1: [3680,3700], 2: [3690,3710], 3: [3700,3720], Index from 1 to 3.
;		default: 1
;
;	wstep - wavelength step size of the output spectra in Angstrom.
;		default: 1.233 angstrom
;
;	pixel - wavelength step size in pixel unit. if pixel keyword is defined, then wstep parameter will be ignored.
;
;	wavearr - User can use user-defined wavelength array for plot. It is used only for plot. It should be regular-spaced.
;
;	ofile - Name of the output file of user spefified data types.
;		default: 'outdata.txt'
;
;
; KEYWORD PARAMETERS: 
;	noplot - if this keyword is setted then there will be no ps file output.
;
;	asahi - ETC also can use theoretical filter throughput from provided by asahi spectra. If this keyword is setted then
;		ftype parameter will be ignored.
;
;	calibration - if this keyword is setted then wavelength calibration light will be used as source. galtemp, inputflux and bbtemp
;		will be ignored. It reads line list from 'calibration_line_list_ori.txt'. 
;
;	calid - User can select calibration source by integer. 0: KrHgNe, 1: Kr, 2: HgNe. Default: KrHgNe
;
;	cspline - if cspline is setted then cspline method will be used for interpolation. Otherwise, linear 
;		interpolation will be used.
;	
;	scflag - if scflag is setted then second order diffracted light will be considered.
;
; OUTPUTS:
;       plot will be exported to 'snr.ps' or other specified file unless noplot keyword is setted
;
; OPTIONAL OUTPUTS:
;	waveout - output wavelength array will be saved in this variable
;
;	snrout - output 'signal to noise' array will be saved in this variable
;
;	signalarr - output signal count will be saved in this variable
;
;	skysignalarr - output sky signal count will be saved in this variable
;
;	noisearr - output noise count will be saved in this variable
;
;	noiseskyarr - output skysignal count will be saved in this variable
;
;	dtypes - Arrays of data type index. User can select data types and export them in ascii file. This parameter should be given
;		in array context, even if only one type will be exported. eg. dtypes=[1],  not dtypes=1
;		user can define output file name by 'ofile' parameter. 
;		index of available data type are given below. (from index no. 1 to 6)
;		1:'wavelength', 2:'snr', 3:'signal', 4:'skysignal', 5:'noise', 6:'noise from sky', 7:'sourceflux'
;
; COMMENTS:
;
; EXAMPLES:
;	1. To check input and output parameters
;	IDL> dotifs_etc, /help
;
;	2. rmag=17, Sc type galaxy at redshift z=0.04, 10 minutes exposure with output name of 'Sc_rmag17_z0.04.ps'
;	IDL> dotifs_etc, exptime=600, band='r', magnitude=17, galtemp='sc', z=0.04, oname='Sc_rmag17_z0.04.ps'
;
;	3. Sc type galaxy with asahi filter, spline interpolation.
;	IDL> dotifs_etc, galtemp='sc', /asahi, /cspline
;
;	4. Sa type galaxy with asahi filter, apply second order contamination, change wavelength bin size as 2.5 pixels.
;	IDL> dotifs_etc, galtemp='sa', /asahi, /scflag, pixel=2.5
;
;	5. Wavelength Calibration source spectrum, change wavelength bin size as 3 angstrom
;	IDL> dotifs_etc, exptime=600, /calibration, calid=0, oname='Calib.ps', wstep=3
;
;	6. Continuum source spectrum
;	IDL> dotifs_etc, exptime=600, bbtemp=5000, oname='BB5000.ps'
;
;	7. Sa type galaxy, plotrange from 5000-6800 angstrom
;	IDL> dotifs_etc, galtemp='sc', plotrange=[5000,6800]
;
;	8. Sa type galaxy, modify plotrange using ptype
;	IDL> dotifs_etc, galtemp='sc', ptype=3
;
;	9. Sa type galaxy, full moon sky
;	IDL> dotifs_etc, galtemp='sc', stype=3
;
;	10. Sa type galaxy, half moon sky, scale sky brightness upto rmag=22
;	IDL> dotifs_etc, galtemp='sc', stype=3, band='r', skymagnitude=22
;	
;	11. Sa type galaxy, extract plotted data
;	IDL> dotifs_etc, galtemp='sc' snrout=snrdata, signalarr=signaldata, noisearr=noisedata, waveout=wavedata
;	IDL> print, wavedata, snrdata, signaldata/noisedata
;
;	12. Sa type galaxy, export some data used for plot. wavelength, signal to noise ratio, and sourceflux. 
;	IDL> dotifs_etc, galtemp='sc', /noplot, ofile='result.txt', dtypes=[1,2,7]
;
;	13. Xenon Arc Lamp Continuum source spectrum
;	IDL> dotifs_etc, exptime=600, /arclamp, oname='Xenon.ps'
;
;
;
;
; MODIFICATION HISTORY:
;	Haeun Chung, 2015, Jun 28, IUCAA, First version
;	Haeun Chung, 2016, Jun 09, IUCAA, Modifed for internal distribution
;	Haeun Chung, 2017, Oct 24, SNU, Add calibration sources. Add Arc lamp source. Fix usd_asahi=1. Add Asachi filter at
;					AOI=10 deg.
;	Haeun Chung, 2018, Nov 12, SNU, Add Littrow ghost. Change gen_cal scheme to add non-zero flux to every pixel
;
function return_littrow_ghost, wave, signal, ccdreflect, cam, g1stR, g1st, g0th
;	lghost=return_littrow_ghost(wave,signal+skysignal, ccdreflect, camtrans, g1stR, g1st, g0th)
	rsub=0.01 ;reflectivity of the grating substrate
	lghost1=total(signal*ccdreflect*cam*cam*g1stR)
	lghost2=total(signal*ccdreflect*cam*cam*g1st*rsub*g0th)
	print, 'Littrow ghost 1 and 2 ', lghost1, lghost2
	lgwave=4800.
	lghost=signal*0

	diff=abs(wave-lgwave)
	sort_idx=sort(diff)
	adjwave=wave[sort_idx[[0,1]]]
	match, adjwave, wave, suba, subb
	binsize=abs(adjwave[0]-adjwave[1])
	tghost=(lghost1+lghost2)
	wave_diff=abs(adjwave[suba]-lgwave)
	wave_ratio=1-wave_diff/binsize
	lghost[subb]=lghost[subb]+tghost*wave_ratio
	return, lghost
end

pro genfilter, params, wave, filter, ftrange, trange
	nwave=n_elements(wave)
	ftrange=double(ftrange)
	trange=double(trange)
	mint=trange[0]
	maxt=trange[1]
	stw=ftrange[0]
	edw=ftrange[1]
	trans=wave*0
	idx=where((wave gt stw) and (wave lt edw), nidx)
	tmpwstep=(edw-stw)/(nidx-1)
	tmpwave=lindgen(nidx)*tmpwstep+stw
	tmpwave2=lindgen(nidx)*tmpwstep+7400d0
	tmptstep=(maxt-mint)/(nidx-1)
	tmpfiltertrans=lindgen(nidx)*tmptstep+mint
	finalfilter=interpol(tmpfiltertrans,tmpwave,wave, spline=cspline)
	finalfilter2=interpol(reverse(tmpfiltertrans),tmpwave2,wave, spline=cspline)
	waveminidx=where(wave le stw)
	wavemaxidx=where(wave ge edw)
	waveidx2=where((wave gt 7400d0) and (wave lt 7400+(edw-stw)))
	waveidx3=where(wave ge 7400+(edw-stw))
	finalfilter[waveminidx]=mint
	finalfilter[wavemaxidx]=maxt
	finalfilter[waveidx2]=finalfilter2[waveidx2]
	finalfilter[waveidx3]=mint
	filter=finalfilter
end


pro dotifs_etc, scflag=scflag, magnitude=magnitude, exptime=exptime, skymagnitude=skymagnitude, plotrange=plotrange, $
	wstep=wstep, ftrange=ftrange, trange=trange, galtemp=galtemp, band=band, z=z, pixel=pixel, help=help, ftype=ftype, ptype=ptype, $
	oname=oname, cspline=cspline, wavearr=wavearr, signalarr=signalarr, noisearr=noisearr, noplot=noplot, inputflux=inputflux, inputwave=inputwave,$
	skysignalarr=skysignalarr, bbtemp=bbtemp, calibration=calibration, wshift=wshift, asahi=asahi, waveout=waveout, noiseskyarr=noiseskyarr, $
	rn=rn, npix_spa=npix_spa, stype=stype, skyfrac=skyfrac, readnoisefrac=readnoisefrac, snrout=snrout, dark=dark, dtypes=dtypes, ofile=ofile, $
	dir=dir, calid=calid, arclamp=arclamp

;system parameters
pri=3.6d0
sec=0.915d0
flag=0
use_asahi=1
ltrghost=1
pixelscale=3700d0/3000   ;angstrom/pixel
dir='~/dotifs_etc/'
if keyword_set(help) then begin
	doc_library, 'etc'
	return
endif
if keyword_set(scflag) then flag=1
if not keyword_set(magnitude) then magnitude=17d0
if not keyword_set(inputparams) then inputparams='input.params'
if not keyword_set(exptime) then exptime=900d0
;if not keyword_set(skymagnitude) then skymagnitude=22d0
if not keyword_set(plotrange) then plotrange=[3700d0,7400d0]
if not keyword_set(wstep) then wstep=3700./3000d0
if not keyword_set(ftrange) then ftrange=[3550d0,3700d0]
if not keyword_set(trange) then trange=[0.001d0, 0.9d0]
if not keyword_set(band) then band='r'
if not keyword_set(ftype) then ftype=1
if keyword_set(pixel) then wstep=pixel*pixelscale
if not keyword_set(pixel) then pixel=wstep/pixelscale
if not keyword_set(rn) then rn=2
if not keyword_set(dark) then dark=0
if not keyword_set(npix_spa) then npix_spa=5
if not keyword_set(stype) then stype=0
if not keyword_set(oname) then filename='snr.ps'
if not keyword_set(dir) then dir='./'
if keyword_set(asahi) then use_asahi=1


rangearr=[$
' (3680-3700, 7400-7420)',$
' (3690-3710, 7400-7420)',$
' (3700-3720, 7400-7420)'$
]


if keyword_set(ftype) then filtername='filter_'+string(ftype, format='(I1)')+rangearr[ftype-1]
if keyword_set(ftype) then begin
case ftype of
	1: ftrange=[3680,3700]
	2: ftrange=[3690,3710]
	3: ftrange=[3700,3720]
endcase
endif

if keyword_set(ptype) then begin
case ptype of
	1: plotrange=[3700,7400]
	2: plotrange=[3700,3760]
	3: plotrange=[7100,7400]
endcase
endif


skysamplingsize=0.4^2*!pi


params=create_struct($
'pridim', 0d,$
'secdim', 0d,$
'objtype', 0,$
'mag', 0d,$
'plotstwave',0d,$
'plotedwave',0d,$
'wavestep',0d,$
'fbstwave',0d,$
'fbedwave',0d,$
'fmaxtrans',0d,$
'fmintrans',0d$
)

;	readparams, inputparams, params		;read telescope parameters
;	params.mag=magnitude
;	pri=params.pridim
;	sec=params.secdim
	stwave=plotrange[0]
	edwave=plotrange[1]
;	wstep=params.wavestep
	telaream2=(pri^2.-sec^2.)/4*!pi		;in m^2
	telarea=telaream2*1e4				;in cm^2
	
	transfile=dir+'trans150626.dat'
	format='D,D,D,D,D,D,D,D,D'
;	readcol, transfile, format=format, wavemicron, skytrans,telmag,col,cam,ccd,g0th,g1st,g2nd, /silent
	readcol, transfile, format=format, wavemicron, skytrans,telmag,col,cam,ccd,g0th_o,g1st_o,g2nd_o, /silent

        file=dir+'vph_160307_new.dat'
        readcol, file, f='D,D', wave_vph_new, eff_vph_new
        eff_vph_new=interpol(eff_vph_new, wave_vph_new/1000, wavemicron)
	g1st=eff_vph_new
	g0th=g1st/g1st_o*g0th_o
	g2nd=g1st/g1st_o*g2nd_o
	g1stR=g1st*5*10.^(-5.)

        file=dir+'altcoating.dat'
        readcol, file, f='D,D,D,D,D,D,D', wave_ac, r0, r20, r30, t0, t20, t30
        eff_caf2_coat=interpol(t0,wave_ac/1000., wavemicron)/100.

        file=dir+'ccd_multi2.txt'
        readcol, file, F=format, wavenm_cm2, ccd_multi2
        waveang_cm2=wavenm_cm2*10
        ccd=interpol(ccd_multi2, waveang_cm2/10000., wavemicron)
	ccdreflect=1-ccd

;        file=dir+'asahi_filter_01_aoi10.dat'
;        readcol, file, F=format, wavenm_filter, asahi_filter_01
;        filter=interpol(asahi_filter_01/100., wavenm_filter/1000., wavemicron)

        col=col/0.995^14.*eff_caf2_coat^14.
        cam=cam/0.995^18.*eff_caf2_coat^18.
	waveang=wavemicron*1d4

if not keyword_set(wavearr) then begin
	nwave=long((waveang[n_elements(waveang)-1]-waveang[0])/wstep)+1
	wave=lindgen(nwave)*wstep+waveang[0]
	diffwave=wave-stwave
	abovezeroidx=where(diffwave ge 0)
	offwave=min(diffwave[abovezeroidx])
	wave=wave-offwave
endif else begin
	nwave=n_elements(wavearr)
	wave=wavearr
	wstep=wavearr[1]-wavearr[0]
	stwave=wave[0]
	edwave=wave[nwave-1]
endelse

	magarr=replicate(magnitude, nwave)

	consth=6.62606957*10^(-34d)
	constc=299792458.

;	photone=!const.h*1.d7*!const.c/(wave*1e-10)
	photone=consth*1.d7*constc/(wave*1e-10)
	sourceflux=mag2flux(magnitude, abwave=wave)
	constflux=replicate(sourceflux[n_elements(sourceflux)/2], nwave)
	sourceflux=constflux

	bandtransfile=dir+band+'filter.dat'
	readcol, bandtransfile, format='D,D', bandwave, bandtrans, /silent

	galtempfile='Flat Magnitude'
if keyword_set(galtemp) then begin
	galtempfile=dir+'kc96/'+galtemp+'_template.ascii'
	readcol, galtempfile, format='D,D', galwave, galflam, /silent

;; apply redshift
	if keyword_set(z) then begin
		galflam=galflam/(1.d0+z)
		galwave=galwave*(1.d0+z)
	endif

	if not keyword_set(z) then z=0.d

	flux2bpmag, bpmag, galflam, galwave, bandtrans, filterwave=bandwave
;print, bpmag
	ratio=10d0^(-0.4d0*(magnitude-bpmag))
	sourceflux=ratio*galflam

	if keyword_set(inputflux) then begin
		if not keyword_set(inputwave) then return
		sourceflux=inputflux
		galwave=inputwave
		galtempfile='User defined input spectrum'
	endif

	sourceflux=interpol(sourceflux, galwave, wave, spline=cspline)
	sourceflux=sourceflux*((wave ge min(galwave)) and (wave le max(galwave)))
endif

;print, constflux
	sourcecount=(sourceflux/photone)*wstep*telarea*exptime*skysamplingsize

if keyword_set(bbtemp) then begin
	sourceflux=planck(wave, bbtemp)*0.015^2.*!pi/(1e5)^2
	sourcecount=planck(wave, bbtemp)*0.015^2.*!pi/(1e5)^2/photone*wstep*exptime
	galtempfile='Blackbody of temperature='+string(bbtemp, format='(F0.1)')+'K'
endif

if keyword_set(arclamp) then begin
	readcol, dir+'Xenon_lamp.txt', format='D,D', xwave, xflam
        sourceflux=interpol(xflam, xwave, wave, spline=cspline)*10.^(-12d)
        sourceflux=sourceflux*((wave ge min(xwave)) and (wave le max(xwave)))
	sourcecount=(sourceflux/photone)*wstep*telarea*exptime*skysamplingsize
	galtempfile='Xenon Arc Lamp'
endif

if keyword_set(calibration) then begin
	calname_arr=['KrHgNe','Kr','HgNe']
if not keyword_set(calid) then calid=0
	calname=calname_arr[calid]
	gen_cal, wave, sourceflux, dir, calname
	sourcecount=sourceflux/photone*wstep*exptime
	galtempfile='Wavelength Calibraiton Source - '+calname
endif



;; sky radiance calculation
;	readcol, 'sky150626.dat', format='D,D', skywave, skyunitcount
;	readcol, dir+'sky150701_newmoon.dat', format='D,D', skywave, skyunitcount, /silent
skyfilearr=['sky150701_newmoon_alt90.dat','sky160311_newmoon_alt45.dat','sky160311_halfmoon.dat','sky160311_fullmoon.dat']
skyfile=skyfilearr[stype]
	readcol, dir+skyfile, format='D,D', skywave, skyunitcount, /silent
	skywave=skywave*10
;	skyphotone=!const.h*1.d7*!const.c/(skywave*1e-10)
	skyphotone=consth*1.d7*constc/(skywave*1e-10)
	skycount=skyunitcount*1d-4*1d-4
	skyflux=skycount*skyphotone
	skyflux=interpol(skyflux, skywave, wave,spline=cspline )
	skyflux=skyflux*((wave ge min(skywave)) and (wave le max(skywave)))

	flux2bpmag, bpskymag, skyflux, wave, bandtrans, filterwave=bandwave

if keyword_set(skymagnitude) then begin
	ratio=10d0^(-0.4d0*(skymagnitude-bpskymag))
	skyflux=ratio*skyflux
endif
	skycount=(skyflux/photone)*wstep*telarea*exptime*skysamplingsize
;pwr=[7350,7450]
;p1st=where((waveang ge pwr[0]) and (waveang le pwr[1]))
;p2nd=where((waveang ge pwr[0]/2) and (waveang le pwr[1]/2))
;print, wave

	ifutrans=0.85d0*0.9   ; 0.9 is just arbitrary factor
;	comtrans=telmag*ifutrans*col*cam*ccd*0.5   ; have no idea why there is 0.5
	comtrans=telmag*ifutrans*col*cam*ccd
	t1st=comtrans*g1st
	t2nd=comtrans*g2nd*0.5   ;light from short wavelength are divided into double wavelength bin
	tsky=skytrans
	t1st=interpol(t1st, waveang, wave, spline=cspline)
	t2nd=interpol(t2nd, waveang, wave, spline=cspline)
	tsky=interpol(tsky, waveang, wave, spline=cspline)

	genfilter, params, wave, filter, ftrange, trange
if use_asahi then begin
	print, 'read asahi filter'
;	readcol, dir+'asahi_filter_trans.dat', format='D,D', afwave, aftrans
	readcol, dir+'asahi_filter_01.dat', format='D,D', afwave, aftrans
;	readcol, dir+'asahi_filter_01_aoi10.dat', format='D,D', afwave, aftrans
	aftrans=aftrans/100.
;	afwave=afwave+wshift-6.6
	afwave=afwave*10.
	filter=interpol(aftrans, afwave, wave, spline=cspline)
	filtername='Asahi filter'
endif
	t1st=t1st*filter
	t2nd=t2nd*filter
;	t2nd=t2nd*filter*0.8

	wave1stidx=where(wave le 7400)
	wave2ndidx=where(wave le 3700)

	wave1stidx=lindgen(nwave)
	wave2ndidx=lindgen(nwave)

;	pwr=[7350,7450]
;	p1st=where((wave ge pwr[0]) and (wave le pwr[1]))
;	p2nd=where((wave ge pwr[0]/2) and (wave le pwr[1]/2))

if (keyword_set(bbtemp)) or (keyword_set(calibration)) or (keyword_set(arclamp)) then begin
	pc1st=t1st*sourcecount
	pc2nd=t2nd*sourcecount
endif else begin
	pc1st=t1st*sourcecount*tsky
	pc2nd=t2nd*sourcecount*tsky
endelse

	skypc1st=t1st*skycount
	skypc2nd=t2nd*skycount

	wave2nd=wave*2
	pc2nd=interpol(pc2nd, wave2nd, wave, spline=cspline)
	skypc2nd=interpol(skypc2nd, wave2nd, wave, spline=cspline)
	pc2nd=pc2nd*((wave ge min(wave2nd)) and (wave le max(wave2nd)))*flag
	skypc2nd=skypc2nd*((wave ge min(wave2nd)) and (wave le max(wave2nd)))*flag

	rn_t=rn*(npix_spa*pixel)^0.5
	dark_t=dark*exptime/3600.*(npix_spa*pixel)^0.5
;print, dark, exptime, npix_spa, pixel, dark_t, rn_t
;	print, rn_t


	signal=pc1st+pc2nd
;	signal=pc1st
	skysignal=skypc1st+skypc2nd

if keyword_set(ltrghost) then begin
	lghost=return_littrow_ghost(wave,signal+skysignal, ccdreflect, cam, g1stR, g1st, g0th)
	signal=signal+lghost
endif

;print, dark, exptime, npix_spa, pixel, dark_t, rn_t, mean((signal+2*skysignal)^0.5)
	noise_poisson=(signal+skysignal+rn_t^2+dark_t^2)^0.5
	noise_sky=(skysignal+rn_t^2+dark_t^2)^0.5
	noise_2nd=pc2nd
	noise=(noise_poisson^2+noise_sky^2)^0.5
;	noise_total=noise_poisson+noise_sky
	nfrac_poisson=signal/noise
	nfrac_sky=skysignal*2/noise
	nfrac_rn=2^0.5*rn_t/noise
	nfrac_dark=2^0.5*dark_t/noise
	snr=signal/noise
	pc2vsntotal=pc2nd/noise
;print, nfrac_rn
;print, rn_t/skysignal^0.5
;print, skysignal^0.5
;print, skysignal[lindgen(500)+100]
;print, n_elements(skysignal)
;print, skycount
	idx=where((wave ge stwave) and (wave le edwave))
	psourceflux=sourceflux[idx]
	pwave=wave[idx]
	psnr=snr[idx]
	psignal=signal[idx]
	pskysignal=skysignal[idx]
	pnoise_sky=noise_sky[idx]
	pnoise=noise[idx]
	pnfrac_poisson=nfrac_poisson[idx]
	pnfrac_sky=nfrac_sky[idx]
	ppc2vsntotal=pc2vsntotal[idx]
	ppc2nd=pc2nd[idx]
	pskyfrac=nfrac_sky[idx]
	prnfrac=nfrac_rn[idx]
	waveout=pwave

	snrout=psnr
	signalarr=psignal
	skysignalarr=pskysignal
	noisearr=pnoise
	noiseskyarr=pnoise_sky
	skyfrac=pskyfrac
	readnoisefrac=prnfrac
	sourcefluxout=psourceflux


if not keyword_set(noplot) then begin	
	plrarr=['','_all','_blueside','_redside']

	if keyword_set(oname) then filename=oname

        !p.thick=3.
        !x.thick=3.
        !y.thick=3.
        !p.charthick=3.
        !p.charsize=1.1

        set_plot, 'ps'
        device, filename=filename, bits_per_pixel=8, /color, xoffset=0, yoffset=1
        device, decomposed=0, xsize=21, ysize=27, /times ;helvetica=1
        xyouts, 0, 0, '!6',/nor


        loadct, 39, /silent
        !p.color=0
        !p.background=255

	charsize=1
	upmar=1.04
	lwmar=0.995
	yposint=0.2
	ysize=yposint
	yst=0.67
	tlen=0.02
	tsty=2
	ttick=1

	i=0
;	xr=[min(pwave)*lwmar, max(pwave)*upmar]
	xr=[min(pwave), max(pwave)]
;	yr=[0, max(psnr)]
	yr=[0, max(psourceflux*1.d17)*upmar]
	pos=[0.15,yst-yposint*i,0.9,yst-yposint*i+ysize]
	gentick, xr, xtickv, xticks, xminor
	gentick, yr, ytickv, yticks, yminor, /exact
        plot, [0,0], [0,0], xrange=xr, yrange=yr, xsty=1, ysty=1, pos=pos, xcharsize=charsize, $
	xticks=xticks, xtickv=xtickv, xminor=xminor, yticks=yticks, ytickv=ytickv, yminor=yminor,$
	/nodata, ytitle=textoidl('f_\lambda (10^{-17} erg/s/cm^2/Ang)'), xtickfor='(A1)',  xticklen=tlen, yticklen=tlen, xgridstyle=tsty, ygridstyle=tsty, /noerase
	for j=0, xticks do begin
		oplot, [xtickv[j],xtickv[j]],[yr[0],yr[1]], linestyle=2, color=0, thick=ttick
	endfor
xgr=[xr,xtickv]
	for j=0, yticks do begin
		oplot, [min(xgr),max(xgr)],[ytickv[j],ytickv[j]], linestyle=2, color=0, thick=ttick
	endfor

	oplot, pwave, psourceflux*1.d17, linestyle=0, color=0

	i=i+1
	xr=[min(pwave), max(pwave)]
	xr=[min(pwave), max(pwave)]
;	yr=[0, max(psnr)]
	yr=[0, max(psnr)*upmar]
	pos=[0.15,yst-yposint*i,0.9,yst-yposint*i+ysize]
	gentick, xr, xtickv, xticks, xminor
	gentick, yr, ytickv, yticks, yminor, /exact
        plot, [0,0], [0,0], xrange=xr, yrange=yr, xsty=1, ysty=1, pos=pos, xcharsize=charsize, $
	xticks=xticks, xtickv=xtickv, xminor=xminor, yticks=yticks, ytickv=ytickv, yminor=yminor,$
	/nodata, ytitle='S/N', xtickfor='(A1)',  xticklen=tlen, yticklen=tlen, xgridstyle=tsty, ygridstyle=tsty, /noerase
	for j=0, xticks do begin
		oplot, [xtickv[j],xtickv[j]],[yr[0],yr[1]], linestyle=2, color=0, thick=ttick
	endfor
xgr=[xr,xtickv]
	for j=0, yticks do begin
		oplot, [min(xgr),max(xgr)],[ytickv[j],ytickv[j]], linestyle=2, color=0, thick=ttick
	endfor
	oplot, pwave, psnr, linestyle=0, color=0
	
	i=i+1
	xr=[min(pwave), max(pwave)]
	yr=[0, max(psignal)*upmar]
;	yr=[0, max(psignal)]
	yint=max(psignal)/5
	pos=[0.15,yst-yposint*i,0.9,yst-yposint*i+ysize]
	gentick, xr, xtickv, xticks, xminor
	gentick, yr, ytickv, yticks, yminor, /exact
        plot, [0,0], [0,0], xrange=xr, yrange=yr, xsty=1, ysty=1, pos=pos, xcharsize=charsize,  $
	xticks=xticks, xtickv=xtickv, xminor=xminor, yticks=yticks, ytickv=ytickv, yminor=yminor,$
	/nodata, ytitle='Signal Count(Black)!CNoise Count(Red)', xtickfor='(A1)',  xticklen=tlen, yticklen=tlen, xgridstyle=tsty, ygridstyle=tsty, /noerase;, xtickint=xint, ytickint=yint
	for j=0, xticks do begin
		oplot, [xtickv[j],xtickv[j]],[yr[0],yr[1]], linestyle=2, color=0,  thick=ttick
	endfor
xgr=[xr,xtickv]
	for j=0, yticks do begin
		oplot, [min(xgr),max(xgr)],[ytickv[j],ytickv[j]], linestyle=2, color=0, thick=ttick
	endfor
	oplot, pwave, psignal, linestyle=0, color=0
	oplot, pwave, pnoise, linestyle=0, color=250


	i=i+1
	xr=[min(pwave), max(pwave)]
	yr=[0, max([max(ppc2nd),max(pnoise)])*upmar]
	yint=0.2
	pos=[0.15,yst-yposint*i,0.9,yst-yposint*i+ysize]
	gentick, xr, xtickv, xticks, xminor
	gentick, yr, ytickv, yticks, yminor, /exact
        plot, [0,0], [0,0], xrange=xr, yrange=yr, xsty=1, ysty=1, pos=pos, xcharsize=charsize, $
	xticks=xticks, xtickv=xtickv, xminor=xminor, yticks=yticks, ytickv=ytickv, yminor=yminor,$
	/nodata, ytitle='2nd Order Signal Count (Blue)!CNoise Count (Red)', xtitle='Wavelength (Angstrom)', $
	 xticklen=tlen, yticklen=tlen, xgridstyle=tsty, ygridstyle=tsty, /noerase
	for j=0, xticks do begin
		oplot, [xtickv[j],xtickv[j]],[yr[0],yr[1]], linestyle=2, color=0, thick=ttick
	endfor
xgr=[xr,xtickv]
	for j=0, yticks do begin
		oplot, [min(xgr),max(xgr)],[ytickv[j],ytickv[j]], linestyle=2, color=0, thick=ttick
;		oplot, [xr[0],xr[1]],[ytickv[j],ytickv[j]], linestyle=2, color=0, thick=ttick
	endfor
;	oplot, pwave, pnfrac_poisson, linestyle=0, color=250
;	oplot, pwave, pnfrac_sky, linestyle=0, color=160
;	oplot, pwave, pnfrac_2nd, linestyle=0, color=50
	oplot, pwave, ppc2nd, linestyle=0, color=50
	oplot, pwave, pnoise, linestyle=0, color=250

;remarks
;if not keyword_set(galtemp) then galtempfile='Flat Magnitude'
if keyword_set(z) then begin
	 zstring=string(z, format='(F8.6)')
endif
if not keyword_set(z) then begin
	 zstring=string(0, format='(F8.6)')
endif
if not keyword_set(galtemp) then zstring=' N/A'

	mag='m_'+band
	remarks=[$
	systime(),$
	'DOTIFS S/N Calculator (ver.09/06/16)',$
	'Galaxy Template: '+galtempfile,$
	'Surface Brightness: '+textoidl(mag)+'='+string(magnitude, format='(I2)')+textoidl('mag/arcsec^2'), $
	'Redshift: z='+zstring, $
	'Exposure: '+string(long(exptime), format='(I4)')+' sec',	$
	textoidl('\Delta\lambda')+':'+string(wstep, format='(F6.3)')+' angstrom ('+string(pixel, format='(F5.2)')+' pixels)'   $
	]	
	nremarks=n_elements(remarks)
	rx=0.9
	ryst=0.973
	ryitv=0.016
	charsize=1
	for i=0, nremarks-1 do begin
		xyouts, rx, ryst-ryitv*i, remarks[i], charsize=charsize, /normal, alignment=1.0
	endfor
	
	xyouts, 0.15, 0.878, filtername, /normal, color=0
;210-10*(ftype-1)

;	export filter information
;	outputfilename='filter_trans_'+string(ftrange[0], format='(I4)')+'_'+string(ftrange[1],format='(I4)')+'.dat'
;	outputfilename='afilter_check.dat'
;	openw, output, outputfilename, /get_lun
;	for i=0l, nwave-1 do begin
;		printf, output, wave[i], filter[i], format='(F11.3,F12.6)'
;	endfor
;	close, output
;	free_lun, output

	device, /close
endif
	qm=string(39B)
	if keyword_set(dtypes) then begin
		if not keyword_set(ofile) then ofile='outdata.txt'
		vararr=['','waveout', 'snrout', 'signalarr', 'skysignalarr', 'noisearr', 'noiseskyarr','sourcefluxout']
		nout=n_elements(dtypes)
		varlist=''
		for j=0, nout-1 do begin
			varlist=varlist+vararr[dtypes[j]]+', '
		endfor
		cmd_string='forprint, '+varlist+'textout='+qm+ofile+qm+', /silent'
	void=execute(cmd_string)
	endif

end
