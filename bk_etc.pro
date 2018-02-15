;+
;pro etc, noflag=noflag, magnitude=magnitude, inputparams=inputparams, exptime=exptime, skymagnitude=skymagnitude, $
;	plotranage=plotrange, wstep=wstep, nwstep=nwstep, ftrange=ftrange, trange=trange, galtemp=galtemp, band=band, $
;	z=z, pixel=pixel, help=help, ftype=ftype, ptype=ptype, cspline=cspline
;-


pro readparams, inputparams, params
	format='A,D,X'
	readcol, inputparams, format=format, parname, value

	idx=where(parname eq 'pridim')
	params.pridim=value[idx]
	idx=where(parname eq 'secdim')
	params.secdim=value[idx]
	idx=where(parname eq 'objtype')
	params.objtype=long(value[idx])
	idx=where(parname eq 'plotstwave')
	params.plotstwave=value[idx]
	idx=where(parname eq 'plotedwave')
	params.plotedwave=value[idx]
	idx=where(parname eq 'wavestep')
	params.wavestep=value[idx]
	idx=where(parname eq 'fbstwave')
	params.fbstwave=value[idx]
	idx=where(parname eq 'fbedwave')
	params.fbedwave=value[idx]
	idx=where(parname eq 'fmaxtrans')
	params.fmaxtrans=value[idx]
	idx=where(parname eq 'fmintrans')
	params.fmintrans=value[idx]

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
;pidx=indgen(35)+215
;pidx=indgen(35)+3240
;print, filter[pidx]
;print, wave[pidx]
end



pro etc, noflag=noflag, magnitude=magnitude, inputparams=inputparams, exptime=exptime, skymagnitude=skymagnitude, plotranage=plotrange, $
	wstep=wstep, nwstep=nwstep, ftrange=ftrange, trange=trange, galtemp=galtemp, band=band, z=z, pixel=pixel, help=help, ftype=ftype, ptype=ptype, $
	oname=oname, cspline=cspline
;ftrange: filter transition region range
;trange: filter transmission min/max
;galtemp: galaxy template name

;system parameters
pri=3.6d0
sec=0.915d0
flag=1
pixelscale=3700d0/3000   ;angstrom/pixel
if keyword_set(help) then begin
	doc_library, 'etc'
	return
endif
if keyword_set(noflag) then flag=0
if not keyword_set(magnitude) then magnitude=17d0
if not keyword_set(inputparams) then inputparams='input.params'
if not keyword_set(exptime) then exptime=900d0
if not keyword_set(skymagnitude) then skymagnitude=22d0
if not keyword_set(plotrange) then plotrange=[3700d0,7400d0]
if not keyword_set(wstep) then wstep=3700./3000d0
if not keyword_set(ftrange) then ftrange=[3550d0,3700d0]
if not keyword_set(trange) then trange=[0.001d0, 0.9d0]
if not keyword_set(band) then band='r'
if keyword_set(pixel) then wstep=pixel*pixelscale
if not keyword_set(pixel) then pixel=wstep/pixelscale

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
	3: plotrange=[7340,7400]
endcase
endif


skysamplingsize=0.4^2*!const.pi


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
	telaream2=(pri^2.-sec^2.)/4*!const.pi		;in m^2
	telarea=telaream2*1e4				;in cm^2
	
	transfile='trans150626.dat'
	format='D,D,D,D,D,D,D,D,D'
	readcol, transfile, format=format, wavemicron, skytrans,telmag,col,cam,ccd,g0th,g1st,g2nd

	waveang=wavemicron*1d4
	nwave=long((waveang[n_elements(waveang)-1]-waveang[0])/wstep)+1
	wave=lindgen(nwave)*wstep+waveang[0]
	diffwave=wave-stwave
	abovezeroidx=where(diffwave ge 0)
	offwave=min(diffwave[abovezeroidx])
	wave=wave-offwave

	magarr=replicate(magnitude, nwave)
	skyarr=replicate(skymagnitude, nwave)

	photone=!const.h*1.d7*!const.c/(wave*1e-10)
	sourceflux=mag2flux(magnitude, abwave=wave)
	constflux=replicate(sourceflux[n_elements(sourceflux)/2], nwave)
	sourceflux=constflux

	bandtransfile=band+'filter.dat'
	readcol, bandtransfile, format='D,D', bandwave, bandtrans
if keyword_set(galtemp) then begin
	galtempfile='./kc96/'+galtemp+'_template.ascii'
	readcol, galtempfile, format='D,D', galwave, galflam
	flux2bpmag, bpmag, galflam, galwave, bandtrans, filterwave=bandwave
	ratio=10d0^(-0.4d0*(magnitude-bpmag))
	sourceflux=ratio*galflam

;	flux2bpmag, bpmag2, sourceflux, galwave, bandtrans, filterwave=bandwave

;print, sourceflux/1d-17
;print, bpmag2

;; apply redshift
if keyword_set(z) then begin
	sourceflux=sourceflux/(1.d0+z)
	galwave=galwave*(1.d0+z)
endif
if not keyword_set(z) then z=0.d
	sourceflux=interpol(sourceflux, galwave, wave, spline=cspline)
	sourceflux=sourceflux*((wave ge min(galwave)) and (wave le max(galwave)))
;print, ratio
;print, 'sourceflux'
;print, sourceflux/1d-17
;print, bpmag
endif

;print, constflux
	sourcecount=(sourceflux/photone)*wstep*telarea*exptime*skysamplingsize


;; sky radiance calculation
;	readcol, 'sky150626.dat', format='D,D', skywave, skyunitcount
	readcol, 'sky150701_newmoon.dat', format='D,D', skywave, skyunitcount
	skywave=skywave*10
	skyphotone=!const.h*1.d7*!const.c/(skywave*1e-10)
	skycount=skyunitcount*1d-4*1d-4
	skyflux=skycount*skyphotone
	skyflux=interpol(skyflux, skywave, wave,spline=cspline )
	skyflux=skyflux*((wave ge min(skywave)) and (wave le max(skywave)))
	skycount=(skyflux/photone)*wstep*telarea*exptime*skysamplingsize

;	skyflux=mag2flux(skyarr, abwave=wave)
;	skycount=(skyflux/photone)*exptime*telarea*wstep*skysamplingsize

;print, n_elements(skyflux)
;print, n_elements(wave)
;print, n_elements(bandtrans)
;print, n_elements(bandwave)
	flux2bpmag, bpskymag, skyflux, wave, bandtrans, filterwave=bandwave
print, bpskymag
;print, total(skycount)
;print, wave
;print, skycount
;idx=lindgen(2000)+1000
;print, skycount[idx]
;print, skycount[idx]
;print, skycount

pwr=[7350,7450]
p1st=where((waveang ge pwr[0]) and (waveang le pwr[1]))
p2nd=where((waveang ge pwr[0]/2) and (waveang le pwr[1]/2))
;print, wave

	ifutrans=0.85d0
;	ifutrans=1
	comtrans=telmag*ifutrans*col*cam*ccd
	t1st=comtrans*g1st
	t2nd=comtrans*g2nd*0.5
	tsky=skytrans
;print, skytrans
;print, t1st[p1st]*0.9
;print, waveang[p1st]

	t1st=interpol(t1st, waveang, wave, spline=cspline)
	t2nd=interpol(t2nd, waveang, wave, spline=cspline)
	tsky=interpol(tsky, waveang, wave, spline=cspline)
;print, t1st
	genfilter, params, wave, filter, ftrange, trange
	t1st=t1st*filter
	t2nd=t2nd*filter

	wave1stidx=where(wave le 7400)
	wave2ndidx=where(wave le 3700)

	wave1stidx=lindgen(nwave)
	wave2ndidx=lindgen(nwave)

;	t1stfilter=wave*0
;	t2ndfilter=wave*0
;	t1stfilter[wave1stidx]=1
;	t2ndfilter[wave2ndidx]=1

;print, t2nd[0:50]
;print, t2ndfilter[0:50]
;	t1st=t1st*t1stfilter
;	t2nd=t2nd*t2ndfilter
;print, t2nd
	pwr=[7350,7450]
	p1st=where((wave ge pwr[0]) and (wave le pwr[1]))
	p2nd=where((wave ge pwr[0]/2) and (wave le pwr[1]/2))
;print, wave[p1st]
;print, ''
;print, t1st[p1st]
;print, ''
;print, t2nd[p2nd]
;print, ''
;print, wave[p2nd]
;print, comtrans[p1st]
;print, filter

;print, n_elements(t2nd), n_elements(sourcecount)
	pc1st=t1st*sourcecount*tsky
	pc2nd=t2nd*sourcecount*tsky
	skypc1st=t1st*skycount
	skypc2nd=t2nd*skycount
;print, ''
;print, sourcecount[p1st]
;print, pc1st[p1st]
;print, ''
;print, sourcecount[p2nd]
;print, pc2nd[p2nd]

;	wave2nd=wave[wave2ndidx]
;	pc2nd=pc2nd[wave2ndidx]
;	skypc2nd=skypc2nd[wave2ndidx]
	wave2nd=wave*2
	pc2nd=interpol(pc2nd, wave2nd, wave, spline=cspline)
	skypc2nd=interpol(skypc2nd, wave2nd, wave, spline=cspline)
	pc2nd=pc2nd*((wave ge min(wave2nd)) and (wave le max(wave2nd)))
	skypc2nd=skypc2nd*((wave ge min(wave2nd)) and (wave le max(wave2nd)))
;print, ''
;print, sourcecount[p1st]
;print, pc1st[p1st]
;print, ''
;print, sourcecount[p2nd]
;print, pc2nd[p1st]
;print, skypc1st
;	skycount=skycount*((wave ge min(stwave)) and (wave le max(edwave)))
	signal=pc1st+pc2nd*flag
	skysignal=skypc1st+skypc2nd
;	signal=pc1st+pc2nd*flag
	noise_poisson=(signal+skysignal)^0.5
	noise_sky=skysignal^0.5
	noise_2nd=pc2nd*flag
;	noise=(noise_poisson^2+noise_sky^2+noise_2nd^2)^0.5
	noise=(noise_poisson^2+noise_sky^2)^0.5
;	noise_total=noise_poisson+noise_sky+noise_2nd
	noise_total=noise_poisson+noise_sky
	nfrac_poisson=noise_poisson/noise_total
	nfrac_sky=noise_sky/noise_total
;	nfrac_2nd=noise_2nd/noise_total
	snr=signal/noise
	pc2vsntotal=pc2nd/noise
;print, noise
;print, skysignal
;print, noise_sky
;print, noise[indgen(200)+200]
;print, wave[indgen(200)+200]
;print, noise_total[indgen(250)+250]
;print, noise_poisson[indgen(250)+250]
;print, noise_2nd[indgen(250)+250]
;print, nfrac_poisson[indgen(250)+250]
;print, snr[p1st]
;print, snr
;print, ''
;print, noise
	idx=where((wave ge stwave) and (wave le edwave))
	psourceflux=sourceflux[idx]
	pwave=wave[idx]
	psnr=snr[idx]
	psignal=signal[idx]
	pnoise=noise[idx]
	pnfrac_poisson=nfrac_poisson[idx]
	pnfrac_sky=nfrac_sky[idx]
;	pnfrac_2nd=nfrac_2nd[idx]
	ppc2vsntotal=pc2vsntotal[idx]
;print, pnfrac_2nd
;print, psnr
	ppc2nd=pc2nd[idx]
;	pnoise2=noise2[idx]
;print, ppc2nd
;print, ''
;print, pnoise2
;print, ''
;print, ppc2nd/pnoise2
;print, pwave
;print, pc2nd[idx]
	
	plrarr=['','_all','_blueside','_redside']

if keyword_set(oname) and keyword_set(ftype) and keyword_set(ptype) then begin
	filename='dotifs_snr_'+galtemp+'_'+band+'mag'+string(magnitude, format='(I2)')+plrarr[ptype]+'_filter'+string(ftype, format='(I1)')+'.ps'
print, filename
endif else begin
	filename='snr.ps'
end
	if flag eq 0 then filename='snr.ps'
        !p.thick=3.
        !x.thick=3.
        !y.thick=3.
        !p.charthick=3.
        !p.charsize=1.1

        set_plot, 'ps'
        device, filename=filename, bits_per_pixel=8, /color, xoffset=0, yoffset=1
        device, decomposed=0, xsize=21, ysize=27, /times ;helvetica=1
        xyouts, 0, 0, '!6',/nor


        loadct, 39
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
	/nodata, ytitle=textoidl('f_\lambda (10^{-17} erg/cm^2/Ang)'), xtickfor='(A1)',  xticklen=tlen, yticklen=tlen, xgridstyle=tsty, ygridstyle=tsty, /noerase
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
if not keyword_set(galtemp) then galtempfile='Flat Magnitude'
if keyword_set(z) then begin
	 zstring=string(z, format='(F8.6)')
endif
if not keyword_set(z) then begin
	 zstring=string(0, format='(F8.6)')
endif
if not keyword_set(galtemp) then zstring=' N/A'

	mag='m_'+band
	remarks=[$
	cpr(),$
	'DOTIFS S/N Calculator',$
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
;	openw, output, outputfilename, /get_lun
;	for i=0l, nwave-1 do begin
;		printf, output, wave[i], filter[i], format='(F11.3,F12.6)'
;	endfor
;	close, output
;	free_lun, output

	device, /close
end
