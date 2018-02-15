;NAME:
;       flux2bpmag, bpmag, flam, flamwave, filtertrans, filterwave=filterwave, errflag=errflag
;
; PURPOSE:
;       Return bandpass AB magnitude
;
; INPUTS: 
;	flam - spectral flux densities per unit wavelength [erg/s/cm2/angstrom]
;	flamwave - wavelength array [Angstrom]
;	filtertrans - filter response expressed as QE. (response per photon)
;
; OPTIONAL INPUTS: 
;	filterwave - in case that filtetrans is given on different wavelength grid of flam, 
;			provide wavelength grid for filter response.
;	errflag - error flag in case of different wavelangth grid.
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS:
;	bpmag - bandpass AB magnitude (double)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;       IDL> flux2bpmag(bpmag, flam, flamwave, filtertrans, filterwave=filterwave, errflag=errflag
;
; MODIFICATION HISTORY:
;	Heun Chung, 2015, Jun 28, IUCAA


;pro flux2bpmag, bpmag, flam, wave, filtertrans_input, filterwave=filterwave, errflag=errflag
pro bk_flux2bpmag, bpmag, flam, wave, filtertrans_input, filterwave=filterwave, errflag=errflag
	mintrans=0.0001
	filtertrans=filtertrans_input
	if keyword_set(filterwave) then begin
		filtertrans=interpol(filtertrans, filterwave, wave)
		filtertrans=(filtertrans ge mintrans)*filtertrans
	endif
	nwave=n_elements(wave)
	if n_elements(filtertrans) ne nwave then begin
		errflag=1
		return
	endif

	filtertrans=double(filtertrans)
	flux=total(flam*filtertrans)
	refflam=replicate(3631.d-23, nwave)/wave^2.*!const.c*1.d10
	refflux=total(refflam*filtertrans)


	bpmag=-2.5*alog10(flux/refflux)
end
