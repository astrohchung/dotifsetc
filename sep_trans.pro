pro sep_trans
	file='trans150626.dat'
	readcol, file, format='D,D,D,D,D,D,D,D,D',micron, sky, telmag, col, cam, ccd, v0th, v1st, v2nd

	names=['sky','telmag','col','cam','ccd','mvph_0th','mvph_1st','mvph_2nd']
	ncomp=n_elements(micron)
	nnames=n_elements(names)
	sarr=dblarr(ncomp, 8)

	sarr[*,0]=sky
	sarr[*,1]=telmag
	sarr[*,2]=col
	sarr[*,3]=cam
	sarr[*,4]=ccd
	sarr[*,5]=v0th
	sarr[*,6]=v1st
	sarr[*,7]=v2nd

	nnames=n_elements(names)
	for i=0, nnames-1 do begin
		oname='./tf/'+names[i]+'.dat'
print, oname
		forprint, micron*10000, reform(sarr[*,i]), textout=oname, /nocomment
	endfor

end
