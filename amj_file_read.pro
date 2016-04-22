;+
; NAME:
;   amj_amj_file_read
;
; PURPOSE:
;   return a piece of sn image from a particular date.
;
; CALLING SEQUENCE:
;   mg = amj_amj_file_read( date [, hdr ] )
;
; INPUTS:
;   date: 	a string in the form yyyy-mm-dd
;   instr: 	Variable that indicates which instrument is used:
;   			1. Reading KPVT-512 files
;   			2. Reading SPMG files
;   			3. Reading MDI files
;   			4. Reading HMI files
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   mg: an mg structure
;
; OPTIONAL OUTPUTS:
;    hdr: the header of the .fits file
;
; MODIFICATION HISTORY:
;
;-

function amj_file_read, datestr, hdr, instr


; Setting file structure-----------------------------------------------------------

; KPVT 512
if instr eq 1 then begin
	dir_str = strmid(datestr,2,2) + strmid(datestr,5,2) + '/'

	fn0 = '/disk/data/munoz/KPVT/' + dir_str
	fn  = fn0 + '*' + strmid(datestr,0,4) + strmid(datestr,5,2) + strmid(datestr,8,2) + '*.fits'
	dn  = file_search( fn, count=cf )
endif

; KPVT SPMG
if instr eq 2 then begin
	dir_str = strmid(datestr,2,2) + strmid(datestr,5,2) + '/'

	fn0 = '/disk/data/munoz/SPMG/' + dir_str
	fn  = fn0 + '*' + strmid(datestr,0,4) + strmid(datestr,5,2) + strmid(datestr,8,2) + '*.fits'
	dn  = file_search( fn, count=cf )
endif		
	
;MDI
if instr eq 3 then begin

	misstol = 0
	avg = 0
	noavg = 0

	dir_str = strmid(datestr,0,4) + '/' + strmid(datestr,5,2) + '/' + strmid(datestr,8,2)

	fn0 = '/disk/data/SOHO/mdi/fd_M_96m_01d/' + dir_str
	fn0 = fn0 + '/fd_M_96m_01d.' + mdi_datestr( datestr )
	
	itry = 0
	IF( keyword_set( ilist ) ) THEN itry = ilist
	m = 0.0
	REPEAT BEGIN
	  stat = 0
	  fn = fn0 + '.' + string( itry, form='(i4.4)' ) + '.fit*'
	  dn = file_search( fn, count=cf )
	  IF( cf GT 0 ) THEN BEGIN
		m = readfits( dn[0], hdr,  /compress )
		mv = sxpar( hdr, 'MISSVALS' )
		ival = sxpar( hdr, 'INTERVAL' )
		IF( mv LE misstol ) THEN stat = 1
		IF( avg AND ( ival LT 150 ) ) THEN stat = 0
		IF( noavg AND ( ival GT 150 ) ) THEN stat = 0
	  ENDIF
	  itry = itry+1
	ENDREP UNTIL( ( stat EQ 1 ) OR ( itry GE 16 ) )	
	
	cf = stat
	
endif


;HMI
if instr eq 4 then begin
	;fn0 = '/Users/santiagovargas/HMI_Magnetograms/'
	fn0 = '/disk/data/munoz/HMI/'
	fn  = fn0 + '*' + strmid(datestr,0,4) + strmid(datestr,5,2) + strmid(datestr,8,2) + '*.fits'
	dn  = file_search( fn, count=cf )


endif

	

IF( cf gt 0 ) THEN begin
	;HMI uses structures for header values
	if instr eq 4 then begin
  
		read_sdo, dn[0], hdr, m, /uncomp_delete
		
		str = hdr.DATE_OBS
		p =  0.0
		b0 = 0.0
		radius = hdr.RSUN_OBS;  use header for radius
		
		x0 = hdr.CRPIX1
        y0 = hdr.CRPIX2

        P0 = hdr.CROTA2


		if ( (P0 ne 0.0) ) then begin
			m =  ROT(m, P0, 1.0, x0, y0, MISSING = !values.f_nan, /cubic)
		endif                
			
	;The rest uses string arrays
	endif else begin
	
		m = readfits( dn[0], hdr)
		str = sxpar( hdr, 'DATE_OBS' )

		; KPVT 512
		if instr eq 1 then begin
			m = m[*,*,2]
			p =  0.0
			b0 = sxpar( hdr, 'EPH_B0' )
			radius = sxpar( hdr, 'EPH_R0' )

			x0 = fxpar(hdr, 'CRPIX1A')
			y0 = fxpar(hdr, 'CRPIX2A')			
		endif

		; KPVT SPMG
		if instr eq 2 then begin
			m = m[*,*,5]
			p =  0.0
			b0 = sxpar( hdr, 'EPH_B0' )
			radius = sxpar( hdr, 'EPH_R0' )/fxpar(hdr,'SCALE')

			x0 = fxpar(hdr, 'CRPIX1A')
			y0 = fxpar(hdr, 'CRPIX2A')			
		endif		
		
		;MDI
		if instr eq 3 then begin

			p =  0.0
			b0 = sxpar( hdr, 'SOLAR_B0' )
			radius = sxpar( hdr, 'OBS_R0' );  use header for radius
			
			x0 = sxpar( hdr, 'X0' ) + sxpar( hdr, 'X_OFFSET' )
			y0 = sxpar( hdr, 'Y0' ) + sxpar( hdr, 'Y_OFFSET' )
			dx = sxpar( hdr, 'XSCALE' )
			dy = sxpar( hdr, 'YSCALE' )
			
			hfx = fxpar(hdr, 'X0');  	Location of the center in x pixels 
			hfy = fxpar(hdr, 'Y0');     Location of the center in y pixels
			P0 = fxpar(hdr, 'P_ANGLE')
			 
			if ( (P0 ne 0.0) ) then begin
			  m =  ROT(m, P0, 1.0, hfx, hfy, MISSING = !values.f_nan, /cubic)
			endif				
		
		endif		
		
	
	endelse
    
    nx =  n_elements( m(*, 0) )
    ny =  n_elements( m(0, *) )
    
    x = float( ( findgen(nx) - x0 ) );  defaults to double!?
    y = float( ( findgen(ny) - y0 ) )
    
    mg = { img:m, x:x, y:y, rad:radius, b0:b0, p:p, date:str}
  	
endif else begin
  mg = 0.0
endelse

return, mg
end
