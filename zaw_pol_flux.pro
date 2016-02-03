pro mdi_pol_mag, mgnt_in, hdr, PF_data, seg_const=seg_const,display=display,prt=prt,iden_save=iden_save,zoom=zoom,hc=hc,ps_file=ps_file,noaa=noaa,file_noaa=file_noaa
;+
; NAME:
;      mdi_faculae
;
; PURPOSE:
;       For a given MDI synoptic white light snapshot, determine the amount of facular
;       pixels
;       
; CATEGORY:1
;       DATA ANALYSIS
;
; CALLING SEQUENCE:
;       mdi_faculae,mgnt_in, hdr, PF_data
;
; INPUTS:
;       mgnt_in: MDI white light
;       hdr: hdr file of the input fits image of the synoptic magnetogram
;
; KEYWORD:
;      prt: display the intages in the different stages of the process
;
; OPTIONAL INPUTS:
;       seg_const: struct holding control parameters for segmentation
;            1. ker_th: kernel pixel threshold
;            2. ar_th: facular region pixel threshold
;            used in region growth
;            3. qr_th: quiet Sun region threshold
;            4. dila_size: dilation structural size, used to emerge
;            the neibouring AR elements into a single AR; default 4 pixels
;            ; used in morphological closing operation
;            5. eros_size: erosion structural size; default 4 pixels; used to remove
;            small AR kernel pixels after the kernel pixel
;            segmentation; used in morphological opending operation
;            6. k_sig (alternate to 1): kernel pixel threshold = k_sig times the background
;standard deviation; default 15
;            7. ar_grow_sig (alternate to 2): active region pixel
;            threshold = ar_grow_sig times background standard deviation; default 5           
;            8. valid_range: value range of valid magnetogram pixels,
;            default=[-20000,20000]; used to deal bad pixels in the iamges
;
; OUTPUTS: 
;       PF_data: A two element array containing the pixels int the north and in the south
; 
;
; OPTIONAL OUTPUT:
;      ptr: print the intermediate information on console
;      iden_save: "features_ars.processes.sav"
;           if iden_save is set, save the images at different processing stages 
;           holding one variable: imgs(7,nx,ny)
;           imgs(0,nx,ny): original input image
;           imgs(1,nx,ny): segmentation of only kernel pixels
;           imgs(2,nx,ny): after removing small and isolated patches
;           using morphological opening operation
;           imgs(3,nx,ny): region growth to find all active region
;           pixels
;           imgs(4,nx,ny): grouping discrete patches to individual
;           active region using dilation,
;           using morphological distance analysis to determine the
;           number of active regions, and the device pixels in each of
;           these active regions
;           imgs(5,nx,ny): showing active region pixels in these
;           identified active regions.
;           imgs(6,nx,ny): showing active region plus dilated pixels
;           after merging active region patches
;
; KEYWORDS:
;            display: display the images in the various processing steps
;            hc:   make hard copy in PS format
;            ps_file: ps file name for showing identification processes
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; Pending improvement:
;   labeling the region number and boundary on the image of segmented
;   active regions
; 
;
; MODIFICATION HISTORY:
;   2011/09/26:  Andres Munoz-Jaramillo:  Adapted from Jie Zhang's AR detection code
;-
;


;define the usage
if n_params() lt 2 then begin
    print,'USAGE:  mdi_segment,mgnt_in,regions'
    return
endif


;Polar Faculae
seg_const_def={ker_th:160.0, ar_th:50.0, qr_th:1.0, dila_size:4.0, eros_size:0.5,k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], gamma:15.0, inte_norm:100.0, pl_cp:70.0, inv_px_tol: 0.95}

if not keyword_set(seg_const) then begin
    seg_const=seg_const_def
endif

;get the input image and dimension
mgnt=mgnt_in
sz=size(mgnt)

;
;set the display window size
;

if keyword_set(display) then begin
    if not keyword_set(zoom) then zoom=0.25
    display_zoom=0.5
    display_xsize=sz[1]*display_zoom & display_ysize=sz[2]*display_zoom  
    
    print_zoom=1.0
    print_xsize=sz[1]*print_zoom & print_ysize=sz[2]*print_zoom  
endif

;Magnetogram image basic information
;
mgnt_valid=where(mgnt ge seg_const.valid_range[0] and mgnt le seg_const.valid_range[1],n_valid,com=ind_invalid,ncom=n_invalid)
if keyword_set(prt) then print,'MAGNETOGRAM: percentage of valid pixels =',string(float(n_valid)/sz[4]*100.0,'(F4.1)')+"%"
mom=moment(mgnt,/nan)
mgnt_mean=mom[0]
if keyword_set(prt) then print,'MAGNETOGRAM: mean =',string(mom[0]^(1.0/1),'(F7.2)')
mgnt_sigma = sqrt(mom[1]) ; the standard deviation of the whole effective image
if keyword_set(prt) then begin
    print,'MAGNETOGRAM: standard deviation =',string((mom[1])^(1.0/2),'(F7.2)')
    print,'MAGNETOGRAM: standdard Skewness =',string((mom[2])^(1.0/3),'(F7.2)')
    print,'MAGNETOGRAM: standard Kurtosis =',string((mom[3])^(1.0/4),'(F7.2)')
endif
;mgnt_med=median(mgnt,/nan)
;if keyword_set(prt) then print,'MAGNETOGRAM: median =',string(im_med,'(F7.2)')
if keyword_set(display) then begin
    display_thresholdu = mgnt_mean+mgnt_sigma*seg_const.k_sig
    display_thresholdd = - display_thresholdu       
 endif

;Define center and radius
hfx = fxpar(hdr, 'X0');35;'CRPIX1');  Location of the center in x pixels 
hfy = fxpar(hdr, 'Y0');+0.3;    Location of the center in y pixels
di = fxpar(hdr,'R_SUN')-6.0;

;Load Solar Coordinates
P0 = fxpar(hdr, 'P_ANGLE')
RD = fxpar(hdr, 'OBS_DIST')/0.0046491
B0 = fxpar(hdr, 'B0')
L0 = fxpar(hdr, 'L0')

;Observer Coordinates
X_scl = fxpar(hdr, 'XSCALE')/60.0
Y_scl = fxpar(hdr, 'YSCALE')/60.0

;Making sure all images have the right north and south

if P0 ne 0.0 then begin

  mgnt =  ROT(mgnt_in, P0, 1.0, hfx, hfy, MISSING = !values.f_nan, /cubic)
  P0 = 0.0    
endif

Xg = transpose(double(floor(findgen(1024,1024)/1024))) - hfx 
Yg = double(floor(findgen(1024,1024)/1024)) - hfy 
R = sqrt( Xg^2 + Yg^2 )

Xg = Xg*X_scl
Yg = Yg*Y_scl

helio = arcmin2hel(Xg[findgen(long(1024.0*1024.0))], Yg[findgen(long(1024.0*1024.0))], p = P0, b0 = B0, r0 = RD, l0 = L0, sphere = 1)

Lath = reform(helio[0,*],1024,1024)
Lonh = reform(helio[1,*],1024,1024)

nps_pl = 100.0
Lon_pl = findgen(nps_pl+1)/nps_pl*360.0 - 180.0

LatN_pl = dblarr(nps_pl+1) + seg_const.pl_cp
LatE_pl = dblarr(nps_pl+1)
LatS_pl = dblarr(nps_pl+1) - seg_const.pl_cp

N_pl = hel2arcmin(LatN_pl, Lon_pl, vsblN, p = P0, b0 = B0, r0 = RD, l0 = L0 , radius = (di+2.0)*X_scl)/X_scl
N_pl_ind = where(vsblN eq 0)
N_pl[0,N_pl_ind] = !values.f_nan
N_pl[1,N_pl_ind] = !values.f_nan


E_pl = hel2arcmin(LatE_pl, Lon_pl, vsblE, p = P0, b0 = B0, r0 = RD, l0 = L0 , radius = (di+2.0)*X_scl)/X_scl
E_pl_ind = where(vsblE eq 0)
E_pl[0,E_pl_ind] = !values.f_nan
E_pl[1,E_pl_ind] = !values.f_nan

S_pl = hel2arcmin(LatS_pl, Lon_pl, vsblS, p = P0, b0 = B0, r0 = RD, l0 = L0 , radius = (di+2.0)*X_scl)/X_scl
S_pl_ind = where(vsblS eq 0)
S_pl[0,S_pl_ind] = !values.f_nan
S_pl[1,S_pl_ind] = !values.f_nan


if keyword_set(display) then begin
    loadct,0,/silent
    window,0,xsize=display_xsize,ysize=display_ysize,retain=2 
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(mgnt,display_xsize,display_ysize),min=display_thresholdd*0.3,max=display_thresholdu*0.3)
    ;loadct, 13
    plots, (N_pl[0,*] + hfx)*display_zoom, (N_pl[1,*] + hfy)*display_zoom ,color='FFFFFF'x,/device,thick=3
    plots, (E_pl[0,*] + hfx)*display_zoom, (E_pl[1,*] + hfy)*display_zoom ,color='FFFFFF'x,/device,thick=3
    plots, (S_pl[0,*] + hfx)*display_zoom, (S_pl[1,*] + hfy)*display_zoom ,color='FFFFFF'x,/device,thick=3
    tvcircle,di*display_zoom, hfx*display_zoom,hfy*display_zoom,color='FFFFFF'x,/device,thick=3
    loadct, 0
endif
if keyword_set(prt) then begin
    set_plot,'PS'
    device, filename='Magnetogram_raw.eps'
    device,/color,bits_per_pixel=8,/portr,/inches,xsize=6.0,ysize=6.0,$
    xoff=0.5,yoff=0.5
    !p.position=[0.0,0.0,1.0,1.0]  
    !x.window=[0.0,1.0]
    !y.window=[0.0,1.0]
    px = !x.window * !d.x_vsize ;Get size of window in device units
    py = !y.window * !d.y_vsize

    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(mgnt,print_xsize,print_ysize),min=display_thresholdd*0.2,max=display_thresholdu*0.2)
    ;loadct, 13
    plots, (N_pl[0,*] + hfx)*print_zoom/sz[1]*px[1], (N_pl[1,*] + hfy)*print_zoom/sz[2]*py[1] ,color='FFFFFF'x,/device,thick=3
    plots, (E_pl[0,*] + hfx)*print_zoom/sz[1]*px[1], (E_pl[1,*] + hfy)*print_zoom/sz[2]*py[1] ,color='FFFFFF'x,/device,thick=3
    plots, (S_pl[0,*] + hfx)*print_zoom/sz[1]*px[1], (S_pl[1,*] + hfy)*print_zoom/sz[2]*py[1] ,color='FFFFFF'x,/device,thick=3
    tvcircle,(di+4.0)*print_zoom/sz[1]*px[1], hfx*print_zoom/sz[1]*px[1],hfy*print_zoom/sz[2]*py[1],color='FFFFFF'x,thick=3
    device, /close
    set_plot,'X'
    loadct, 0, /silent           
endif
imgs0=mgnt

;Correct line of sight magnetic field assuming field is radial

;Observers' unit vector
Xobs = cos(B0*!dtor)*cos(L0*!dtor)
Yobs = cos(B0*!dtor)*sin(L0*!dtor)
Zobs = sin(B0*!dtor)

M_corr = cos(Lath*!dtor)*cos(Lonh*!dtor)*Xobs + cos(Lath*!dtor)*sin(Lonh*!dtor)*Yobs + sin(Lath*!dtor)*Zobs
mgnt_crr = mgnt/M_corr

if keyword_set(display) then begin
    loadct,0,/silent
    if not keyword_set(ps) then window,4,xsize=display_xsize,ysize=display_ysize,retain=2 
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(mgnt_crr,display_xsize,display_ysize),min=display_thresholdd*0.7,max=display_thresholdu*0.7)
    loadct, 13
    plots, (N_pl[0,*] + hfx)*display_zoom, (N_pl[1,*] + hfy)*display_zoom ,color=200/6,/device,thick=3
    plots, (E_pl[0,*] + hfx)*display_zoom, (E_pl[1,*] + hfy)*display_zoom ,color=200*3/6,/device,thick=3
    plots, (S_pl[0,*] + hfx)*display_zoom, (S_pl[1,*] + hfy)*display_zoom ,color=200/6,/device,thick=3
    tvcircle,di*display_zoom, hfx*display_zoom,hfy*display_zoom,color=200*2/6,/device,thick=3
    loadct, 0
endif
imgs1=mgnt_crr


;Calculate the Element of Area

;Converting to Cartesian
Xar = cos(Lath*!dtor)*cos(Lonh*!dtor)
Yar = cos(Lath*!dtor)*sin(Lonh*!dtor)
Zar = sin(Lath*!dtor)

;Finding the first diagonal vector
Parx = dblarr(sz[1],sz[2])
Parx[1:sz[1]-2,1:sz[2]-2] = Xar[0:sz[1]-3,0:sz[2]-3] - Xar[2:sz[1]-1,2:sz[2]-1]
Pary = dblarr(sz[1],sz[2])
Pary[1:sz[1]-2,1:sz[2]-2] = Yar[0:sz[1]-3,0:sz[2]-3] - Yar[2:sz[1]-1,2:sz[2]-1]
Parz = dblarr(sz[1],sz[2])
Parz[1:sz[1]-2,1:sz[2]-2] = Zar[0:sz[1]-3,0:sz[2]-3] - Zar[2:sz[1]-1,2:sz[2]-1]

P = [reform(Parx,1,sz[1]*sz[2]), reform(Pary,1,sz[1]*sz[2]), reform(Parz,1,sz[1]*sz[2])]

;Finding the second diagonal vector
Qarx = dblarr(sz[1],sz[2])
Qarx[1:sz[1]-2,1:sz[2]-2] = Xar[2:sz[1]-1,0:sz[2]-3] - Xar[0:sz[1]-3,2:sz[2]-1]
Qary = dblarr(sz[1],sz[2])
Qary[1:sz[1]-2,1:sz[2]-2] = Yar[2:sz[1]-1,0:sz[2]-3] - Yar[0:sz[1]-3,2:sz[2]-1]
Qarz = dblarr(sz[1],sz[2])
Qarz[1:sz[1]-2,1:sz[2]-2] = Zar[2:sz[1]-1,0:sz[2]-3] - Zar[0:sz[1]-3,2:sz[2]-1]

Q = [reform(Qarx,1,sz[1]*sz[2]), reform(Qary,1,sz[1]*sz[2]), reform(Qarz,1,sz[1]*sz[2])]

Atmp = crosspn(P,Q)

mgnt_ar = reform(sqrt(Atmp[0,*]^2.0 + Atmp[1,*]^2.0 + Atmp[2,*]^2.0)*6.955e10*6.955e10/8.0, sz[1], sz[2])
ind_npf = where(R gt di)
mgnt_ar(ind_npf) = 0.0

if keyword_set(display) then begin
    loadct,0,/silent
    if not keyword_set(ps) then window,2,xsize=display_xsize,ysize=display_ysize,retain=2 
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(mgnt_ar,display_xsize,display_ysize),min=min(mgnt_ar,/nan),max=max(mgnt_ar,/nan))
    loadct, 13
    plots, (N_pl[0,*] + hfx)*display_zoom, (N_pl[1,*] + hfy)*display_zoom ,color=200,/device,thick=3
    plots, (E_pl[0,*] + hfx)*display_zoom, (E_pl[1,*] + hfy)*display_zoom ,color=200*3/6,/device,thick=3
    plots, (S_pl[0,*] + hfx)*display_zoom, (S_pl[1,*] + hfy)*display_zoom ,color=200,/device,thick=3
    tvcircle,di*display_zoom, hfx*display_zoom,hfy*display_zoom,color=200*2/6,/device,thick=3
    loadct, 0
endif
imgs2=mgnt_ar

;Calculating flux

mgnt_flx = mgnt_ar*mgnt_crr

if keyword_set(display) then begin
    loadct,0,/silent
    if not keyword_set(ps) then window,6,xsize=display_xsize,ysize=display_ysize,retain=2 
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(mgnt_flx,display_xsize,display_ysize),min=min(mgnt_flx/10.0,/nan),max=max(mgnt_flx/10.0,/nan))
    loadct, 13
    plots, (N_pl[0,*] + hfx)*display_zoom, (N_pl[1,*] + hfy)*display_zoom ,color=200,/device,thick=3
    plots, (E_pl[0,*] + hfx)*display_zoom, (E_pl[1,*] + hfy)*display_zoom ,color=200*3/6,/device,thick=3
    plots, (S_pl[0,*] + hfx)*display_zoom, (S_pl[1,*] + hfy)*display_zoom ,color=200,/device,thick=3
    tvcircle,di*display_zoom, hfx*display_zoom,hfy*display_zoom,color=200*2/6,/device,thick=3
    loadct, 0
endif
imgs3=mgnt_flx


;North
n_swt = 0.0
pc_pxindn = where(R le di and Lath ge seg_const.pl_cp, npc_pxn);   Polar cap pixels
vpc_pxindn = where(R le di and finite(mgnt) and Lath ge seg_const.pl_cp, vnpc_pxn);  Valid pixels inside the polar cap
pc_posn = where(R le di and Lath ge seg_const.pl_cp and mgnt_crr gt 0.0, npc_posn); Pixels with possitive corrected field inside the polar cap
pc_negn = where(R le di and Lath ge seg_const.pl_cp and mgnt_crr lt 0.0, npc_negn); Pixels with negative corrected field inside the polar cap

if ( npc_negn eq 0 or npc_posn eq 0 ) then begin
    print, 'Northern hemisphere no polarity mixture' 
    n_swt = 3.0
endif    
if ( float(vnpc_pxn)/float(npc_pxn) le seg_const.inv_px_tol or vpc_pxindn[0] eq -1 ) then begin
    print, 'Northen hemisphere has more than ' + string((1.0-seg_const.inv_px_tol)*100.0) + '% invalid pixels' 
    n_swt = 2.0
endif
if ( n_swt ne 0.0 ) then begin 

    intf_n = !values.f_nan ;Uncorrected Mean polar cap field
    intfc_n = !values.f_nan ;Corrected Mean polar cap field
    unsf_n = !values.f_nan ;Total polar unsigned flux
    sf_n = !values.f_nan ;Total polar signed flux
    posf_n = !values.f_nan ;Total positive polar flux
    negf_n = !values.f_nan ;Total negative polar flux
    vnpc_pxn = !values.f_nan
    visa_n = !values.f_nan
        
endif else begin   
        
    intf_n = mean( mgnt(vpc_pxindn),/double, /nan) ;Uncorrected Mean polar cap field
    intfc_n = mean( mgnt_crr(vpc_pxindn),/double, /nan) ;Corrected Mean polar cap field
    unsf_n = total( abs(mgnt_flx(vpc_pxindn)),/double, /nan) ;Total polar unsigned flux
    sf_n = total( mgnt_flx(vpc_pxindn),/double, /nan) ;Total polar signed flux
    posf_n = total( mgnt_flx(pc_posn),/double, /nan) ;Total positive polar flux
    negf_n = total( mgnt_flx(pc_negn),/double, /nan) ;Total negative polar flux
    visa_n = total( mgnt_ar(vpc_pxindn), /double, /nan) ;Total visible area
            
endelse


;South
s_swt = 0.0
pc_pxinds = where(R le di and Lath le -seg_const.pl_cp, npc_pxs);   Polar cap pixels
vpc_pxinds = where(R le di and finite(mgnt) and Lath le -seg_const.pl_cp, vnpc_pxs);  Valid pixels inside the polar cap
pc_poss = where(R le di and Lath le -seg_const.pl_cp and mgnt_crr gt 0.0, npc_poss); Pixels with possitive corrected field inside the polar cap
pc_negs = where(R le di and Lath le -seg_const.pl_cp and mgnt_crr lt 0.0, npc_negs); Pixels with negative corrected field inside the polar cap

if ( npc_negs eq 0 or npc_poss eq 0 ) then begin
    print, 'Southern hemisphere no polarity mixture' 
    s_swt = 3.0
endif
if ( float(vnpc_pxs)/float(npc_pxs) le seg_const.inv_px_tol or vpc_pxinds[0] eq -1 ) then begin
    print, 'Southern hemisphere has more than ' + string((1.0-seg_const.inv_px_tol)*100.0) + '% invalid pixels' 
    s_swt = 2.0
endif

    
if ( s_swt ne 0.0 ) then begin 

    intf_s = !values.f_nan ;Uncorrected Mean polar cap field
    intfc_s = !values.f_nan ;Corrected Mean polar cap field
    unsf_s = !values.f_nan ;Total polar unsigned flux
    sf_s = !values.f_nan ;Total polar signed flux
    posf_s = !values.f_nan ;Total positive polar flux
    negf_s = !values.f_nan ;Total negative polar flux
    vnpc_pxs = !values.f_nan
    visa_s = !values.f_nan
    
endif else begin   
        
    intf_s = mean( mgnt(vpc_pxinds),/double, /nan) ;Uncorrected Mean polar cap field
    intfc_s = mean( mgnt_crr(vpc_pxinds),/double, /nan) ;Corrected Mean polar cap field
    unsf_s = total( abs(mgnt_flx(vpc_pxinds)),/double, /nan) ;Total polar unsigned flux
    sf_s = total( mgnt_flx(vpc_pxinds),/double, /nan) ;Total polar signed flux
    posf_s = total( mgnt_flx(pc_poss),/double, /nan) ;Total positive polar flux
    negf_s = total( mgnt_flx(pc_negs),/double, /nan) ;Total negative polar flux
    visa_s = total( mgnt_ar(vpc_pxinds), /double, /nan) ;Total visible area
            
endelse

PF_data = [intf_n, intfc_n, unsf_n, sf_n, posf_n, negf_n, vnpc_pxn, visa_n, n_swt, intf_s, intfc_s, unsf_s, sf_s, posf_s, negf_s, vnpc_pxs, visa_s, s_swt] 
 
return
end
