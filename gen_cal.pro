pro gen_cal, wave_in, flux_out, dir, calname
;	readcol, dir+'calibration_line_list_ori.txt', format='D,D,X', lwave, lflux
	readcol, dir+'calibration_line_list_'+calname+'_171024.txt', format='D,D,X', lwave, lflux
	wave=wave_in


;	off=1
;	lwave_all=[]
;	for i=0, 4 do begin
;		lwave_all=[lwave_all,(lwave+off*(i-2))]
;	endfor

;print, lwave
;print, lflux
;	pwave=(lindgen(50000l)+30000d)/10
	

	nline=n_elements(lwave)

	flux_out=wave*0.d
	for i=0, nline-1 do begin
		if (lwave[i] le wave[0]) or (lwave[i] ge wave[n_elements(wave)-1]) then continue
		diff=abs(wave-lwave[i])
		sort_idx=sort(diff)
		adjwave=wave[sort_idx[[0,1]]]
		match, adjwave, wave, suba, subb
		binsize=abs(adjwave[0]-adjwave[1])
		nflux=lflux[i]/binsize   ; to make flux per ang
;		wave_diff=abs(adjwave-lwave[i])
		wave_diff=abs(adjwave[suba]-lwave[i])
		wave_ratio=1-wave_diff/binsize
		flux_out[subb]=flux_out[subb]+nflux*wave_ratio*1.d-8
;		flux_out[subb]=nflux*wave_ratio*1.d-8
;		flux_out[subb]=nflux*wave_ratio*1.d-10
	endfor

	zero_idx=where(0.d eq flux_out, zerocount)
;print, zerocount
;	flux_out[zero_idx]=dblarr(zerocount)+nflux*wave_ratio[0]*1.d-8*0.01
	flux_out=flux_out+nflux*wave_ratio[0]*1.d-8*0.01

;print, flux_out[lindgen(50)]
;print, wave_in
;stop
end
