pro zaw_pol_flux, start_date, end_date, instr, seg_const=seg_const, display=display
; NAME
;      zaw_pol_flux
;
; PURPOSE:
;       For a given MDI synoptic white light snapshot, determine the amount of facular
;       pixels
;       
; CATEGORY:1
;       DATA ANALYSIS
;
; CALLING SEQUENCE:
;       zaw_pol_flux,
;
; INPUTS:
;       
;
; KEYWORD:
;      prt: display the intages in the different stages of the process
;
; OPTIONAL INPUTS:
;       seg_const: struct holding control parameters for segmentation
;            1. ker_th: kernel pixel threshold
;            2. ar_th: facular region pixel threshold
;                   used in region growth
;            3. qr_th: quiet Sun region threshold
;            4. dila_size: dilation structural size, used to emerge
;                   the neibouring AR elements into a single AR; default 4 pixels
;                   used in morphological closing operation
;            5. eros_size: erosion structural size; default 4 pixels; used to remove
;                   small AR kernel pixels after the kernel pixel
;                   segmentation; used in morphological opending operation
;            6. k_sig (alternate to 1): kernel pixel threshold = k_sig times the background
;                   standard deviation; default 15
;            7. ar_grow_sig (alternate to 2): active region pixel
;                   threshold = ar_grow_sig times background standard deviation; default 5           
;            8. valid_range: value range of valid magnetogram pixels,
;                   default=[-20000,20000]; used to deal bad pixels in the iamges
;   
; OUTPUTS: 
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
;   2011/09/26: Andres Munoz-Jaramillo:  Adapted from Jie Zhang's AR detection code
;   2016/02/06: Zach Werginz: Adapted to work with amj_coord.pro and amj_file_read.pro
;

;Define the usage
if n_params() lt 2 then begin
    print,'USAGE:  zaw_pol_flux, start_date, end_date, instr'
    print, "instr: 1 = KPVT 512"
    print, "instr: 2 = KPVT SPMG"
    print, "instr: 3 = MDI"
    print, "instr: 4 = HMI"
    return
endif

deg_lim = 65.0
inv_px_tol = 0.85
stat = 1

;Set instrument parameters
;First reference day for keeping track time easily
if instr eq 1 then DayOff = julday(1,1,1970);   KPVT 512
if instr eq 2 then DayOff = julday(1,1,1970);   KPVT SPMG
if instr eq 3 then DayOff = julday(1,1,1993);   MDI
if instr eq 4 then DayOff = julday(1,1,2009);   HMI

mdi_i  = julday(strmid(start_date,5,2),strmid(start_date,8,2),strmid(start_date,0,4))-DayOff
mdi_f  = julday(strmid(end_date,5,2),strmid(end_date,8,2),strmid(end_date,0,4))-DayOff
PF_data = {pf_day, mdi_i: 0L, date: '', intf_n: !values.f_nan, intfc_n: !values.f_nan, unsf_n: !values.f_nan, sf_n: !values.f_nan, posf_n: !values.f_nan, $
            negf_n: !values.f_nan, vnpc_pxn: !values.f_nan, visa_n: !values.f_nan, max_n: !values.f_nan, n_swt: !values.f_nan, intf_s: !values.f_nan, intfc_s: !values.f_nan,$
            unsf_s: !values.f_nan, sf_s: !values.f_nan, posf_s: !values.f_nan, negf_s: !values.f_nan, vnpc_pxs: !values.f_nan, visa_s: !values.f_nan, max_s: !values.f_nan, s_swt: !values.f_nan}

REPEAT begin
    
    ;Reading files 
    repeat begin
            caldat, mdi_i + DayOff, Month, Day, Year
            date = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)
            print, date
            mg = amj_file_read( date, hdr, instr )
        
            s = size(mg)    
            
            mdi_i = mdi_i + 1          
    endrep until (s[2] eq 8)

    if (date eq end_date) then stat = 0

    ;Reading header variables

    ;HMI uses structures for header values
    if instr eq 4 then begin
        date = hdr.DATE_OBS

        ;Define center and radius
        hfx = hdr.CRPIX1 ;  Location of the center in x pixels 
        hfy = hdr.CRPIX2 ;    Location of the center in y pixels
        di = hdr.RSUN_OBS/hdr.CDELT1;

        ;Load Solar Coordinates
        P0 = 0.0
        RD = hdr.DSUN_OBS/hdr.RSUN_REF
        B0 = hdr.CRLT_OBS
        L0 = hdr.CRLN_OBS

        ;Observer Coordinates
        X_scl = hdr.CDELT1/60.0
        Y_scl = hdr.CDELT2/60.0
        
    endif else begin


        date = fxpar(hdr, 'DATE_OBS')

        ;KPVT-512
        if instr eq 1 then begin
        
            ;Define center and radius
            hfx = fxpar(hdr, 'CRPIX1A');35;'CRPIX1');  Location of the center in x pixels 
            hfy = fxpar(hdr, 'CRPIX2A');+1.0;    Location of the center in y pixels
            di = fxpar(hdr,'EPH_R0');

            ;Load Solar Coordinates
            P0 = 0.0
            RD = !values.f_nan
            B0 = fxpar(hdr, 'EPH_B0')
            L0 = fxpar(hdr, 'EPH_L0')

            ;Observer Coordinates
            X_scl = fxpar(hdr, 'CDELT1')*fxpar(hdr, 'CRR_SCLX')/60.0
            Y_scl = fxpar(hdr, 'CDELT2')*fxpar(hdr, 'CRR_SCLY')/60.0

        endif

        ;MDI
        if instr eq 3 then begin
        
            ;Define center and radius
            hfx = fxpar(hdr, 'X0');  Location of the center in x pixels 
            hfy = fxpar(hdr, 'Y0');  Location of the center in y pixels
            di = fxpar(hdr,'R_SUN');

            ;Load Solar Coordinates
            P0 = 0.0
            RD = fxpar(hdr, 'OBS_DIST')/0.0046491
            B0 = fxpar(hdr, 'B0')
            L0 = fxpar(hdr, 'L0')

            ;Observer Coordinates
            X_scl = fxpar(hdr, 'XSCALE')/60.0
            Y_scl = fxpar(hdr, 'YSCALE')/60.0    
        
        endif

    endelse

    ;Calculate heliospheric data    
    amj_coord, mg.img, hdr, CRD, instr
    s = size(CRD.im_raw)
    Xg = transpose(double(floor(findgen(s[1],s[1])/s[1]))) - hfx 
    Yg = double(floor(findgen(s[1],s[1])/s[1])) - hfy 
    R = sqrt( Xg^2 + Yg^2 )*X_scl*60.0
    
    ;;North
    n_swt = 0.0
    pc_pxindn = where(R le di and CRD.Lath ge deg_lim, npc_pxn);   Polar cap pixels
    vpc_pxindn = where(R le di and finite(CRD.im_raw) and CRD.Lath ge deg_lim, vnpc_pxn);  Valid pixels inside the polar cap
    pc_posn = where(R le di and CRD.Lath ge deg_lim and CRD.im_crr gt 0.0, npc_posn); Pixels with positive corrected field inside the polar cap
    pc_negn = where(R le di and CRD.Lath ge deg_lim and CRD.im_crr lt 0.0, npc_negn); Pixels with negative corrected field inside the polar cap

    if ( npc_negn eq 0 or npc_posn eq 0 ) then begin    
        print, 'Northern hemisphere no polarity mixture' 
        n_swt = 3.0
    endif    
    if ( float(vnpc_pxn)/float(npc_pxn) le inv_px_tol or vpc_pxindn[0] eq -1 ) then begin
        print, 'Northen hemisphere has more than ' + string((1.0-inv_px_tol)*100.0) + '% invalid pixels' 
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
        max_n = !values.f_nan ;Maximum pixel flux            
    endif else begin   
            
        intf_n = mean( CRD.im_raw(vpc_pxindn),/double, /nan) ;Uncorrected Mean polar cap field
        intfc_n = mean( CRD.im_crr(vpc_pxindn),/double, /nan) ;Corrected Mean polar cap field
        unsf_n = total( abs(CRD.mgnt_flx(vpc_pxindn)),/double, /nan) ;Total polar unsigned flux
        sf_n = total( CRD.mgnt_flx(vpc_pxindn),/double, /nan) ;Total polar signed flux
        posf_n = total( CRD.mgnt_flx(pc_posn),/double, /nan) ;Total positive polar flux
        negf_n = total( CRD.mgnt_flx(pc_negn),/double, /nan) ;Total negative polar flux
        visa_n = total( CRD.mgnt_ar(vpc_pxindn), /double, /nan) ;Total visible area
        max_n = max(abs( CRD.mgnt_flx(vpc_pxindn) ), /nan) ;Maximum pixel flux 
                
    endelse


    ;South
    s_swt = 0.0
    pc_pxinds = where(R le di and CRD.Lath le -deg_lim, npc_pxs);   Polar cap pixels
    vpc_pxinds = where(R le di and finite(CRD.im_raw) and CRD.Lath le -deg_lim, vnpc_pxs);  Valid pixels inside the polar cap
    pc_poss = where(R le di and CRD.Lath le -deg_lim and CRD.im_crr gt 0.0, npc_poss); Pixels with possitive corrected field inside the polar cap
    pc_negs = where(R le di and CRD.Lath le -deg_lim and CRD.im_crr lt 0.0, npc_negs); Pixels with negative corrected field inside the polar cap

    if ( npc_negs eq 0 or npc_poss eq 0 ) then begin
        print, 'Southern hemisphere no polarity mixture' 
        s_swt = 3.0
    endif
    if ( float(vnpc_pxs)/float(npc_pxs) le inv_px_tol or vpc_pxinds[0] eq -1 ) then begin
        print, 'Southern hemisphere has more than ' + string((1.0-inv_px_tol)*100.0) + '% invalid pixels' 
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
        max_s = !values.f_nan ;Maximum pixel flux 
        
    endif else begin   
            
        intf_s = mean( CRD.im_raw(vpc_pxinds),/double, /nan) ;Uncorrected Mean polar cap field
        intfc_s = mean( CRD.im_crr(vpc_pxinds),/double, /nan) ;Corrected Mean polar cap field
        unsf_s = total( abs(CRD.mgnt_flx(vpc_pxinds)),/double, /nan) ;Total polar unsigned flux
        sf_s = total( CRD.mgnt_flx(vpc_pxinds),/double, /nan) ;Total polar signed flux
        posf_s = total( CRD.mgnt_flx(pc_poss),/double, /nan) ;Total positive polar flux
        negf_s = total( CRD.mgnt_flx(pc_negs),/double, /nan) ;Total negative polar flux
        visa_s = total( CRD.mgnt_ar(vpc_pxinds), /double, /nan) ;Total visible area
        max_s = max(abs( CRD.mgnt_flx(vpc_pxinds) ), /nan) ;Maximum pixel flux 
                
    endelse
    tmp_PF = {pf_day, mdi_i, date, intf_n, intfc_n, unsf_n, sf_n, posf_n, negf_n, vnpc_pxn, visa_n, max_n, n_swt, intf_s, intfc_s, unsf_s, sf_s, posf_s, negf_s, vnpc_pxs, visa_s, max_s, s_swt}
    PF_data = [PF_data, tmp_PF]
ENDREP UNTIL (mdi_i gt mdi_f)

filename = string('PF_data' + start_date + '_' + end_date + '.csv')
write_csv, filename, PF_data

END
