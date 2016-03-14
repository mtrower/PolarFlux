pro zaw_pol_flux_batch, start_date, end_date, instr
; NAME:
;	   zaw_pol_flux_batch
;
; PURPOSE:  To calculate the polar flux bounded by a specific latitude
;               in a batch process through several magnetograms
;   
; CALLING SEQUENCE: zaw_pol_flux_batch
;       
; INPUTS:  none 
;   
; OPTIONAL INPUTS:
;   
; KEYWORDS:
;
; MODIFICATION HISTORY:
;

PF_data = {pf_day, mdi_i: 0L, date: '', intf_n: !values.f_nan, intfc_n: !values.f_nan, unsf_n: !values.f_nan, sf_n: !values.f_nan, posf_n: !values.f_nan, $
            negf_n: !values.f_nan, vnpc_pxn: !values.f_nan, visa_n: !values.f_nan, max_n: !values.f_nan, n_swt: !values.f_nan, intf_s: !values.f_nan, intfc_s: !values.f_nan,$
            unsf_s: !values.f_nan, sf_s: !values.f_nan, posf_s: !values.f_nan, negf_s: !values.f_nan, vnpc_pxs: !values.f_nan, visa_s: !values.f_nan, max_s: !values.f_nan, s_swt: !values.f_nan}
inv_px_tol = 1.00

for i = 0, 100 do begin
	zaw_pol_flux, start_date, end_date, instr, PF, inv_px_tol = inv_px_tol - i*.01
	tmp_PF = {pf_day, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
    PF_data = [PF_data, tmp_PF]
endfor

End