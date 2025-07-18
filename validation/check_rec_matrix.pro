dir = ROUTINE_DIR()
params = (read_params_file(dir+'params_alasca_leo_sh.pro',/expand))[0]

rec = readfits('/home/bstadler/passata/ALASCA_LEO/calibration_sh_Tenerife_elevation15deg_PAA25murad_/rec/ALASCA_LEO_ps120p0.008_shs10x10_wl589_fv10.0_np24_th0.25_mn5_ce.fits', ext=1)
p8 = diag_matrix(matrix_multiply(transpose(rec),(rec)))

rec = readfits('/home/bstadler/passata/ALASCA_LEO/calibration_sh_Tenerife_elevation15deg_PAA25murad/rec/ALASCA_LEO_ps120p0.008_shs10x10_wl589_fv10.0_np24_th0.25_mn8_ce.fits', ext=1)
p9 = diag_matrix(matrix_multiply(transpose(rec),(rec)))

;rec = readfits('/home/bstadler/passata/ALASCA_LEO/calibration_sh/rec/ALASCA_LEO_ps120p0.008_shs10x10_wl589_fv10.0_np24_th0.25_mn8_ce.fits', ext=1)
;p19 = diag_matrix(matrix_multiply(transpose(rec),(rec)))

;rec = readfits('/home/bstadler/passata/ALASCA_LEO/calibration_pyr/rec/ALASCA_LEO_ps120p0.008_pyr40x40_wl589_fv5.0_ft3.0_ma0_bn1_th0.30a0.30b_mn8.fits', ext=1)
;pyr8 = diag_matrix(matrix_multiply(transpose(rec),(rec)))
;
;rec = readfits('/home/bstadler/passata/ALASCA_LEO/calibration_pyr/rec/ALASCA_LEO_ps120p0.008_pyr40x40_wl589_fv5.0_ft3.0_ma0_bn1_th0.30a0.30b_mn19.fits', ext=1)
;pyr19 = diag_matrix(matrix_multiply(transpose(rec),(rec)))

;im_sh = readfits('/home/bstadler/passata/ALASCA_LEO/calibration_pyr_Tenerife_elevation30deg_PAA50murad/im/ALASCA_LEO_ps120p0.008_shs12x12_wl589_fv10.0_np20_th0.25_mn66_ce.fits', ext=1)
;s = (size(im_sh,/dim))[1]
scale_sh = 1 ;abs(mean(reform(im_sh[1,0:s/2.-1])))

;im_pyr = readfits('/home/bstadler/passata/ALASCA_LEO/calibration_pyr/im/ALASCA_LEO_ps120p0.008_pyr40x40_wl589_fv5.0_ft3.0_ma0_bn1_th0.30a0.30b_mn19.fits', ext=1)
;s = (size(im_pyr,/dim))[1]
;scale_pyr = mean(reform(im_pyr[1,0:s/2.-1]))
;SVDC, rec, W, U, V

;window, 7, xs=1000, ys=600
vcol = [1l,255l]
p = plot(sqrt(p8), "b-2", XRANGE=[0,16], YRANGE=[0,max(sqrt(p8))], TITLE='SH WFS', YTITLE='noise propagation coefficient [nm]', XTITLE='mode number', dim=[1000,600], NAME='20 cm')
p2 = plot(sqrt(p9), "r-2", XRANGE=[0,16], /OVERPLOT, NAME='30 cm')
leg = legend(TARGET=[p,p2],POSITION=[0.9,0.8])

;p = plot(sqrt(pyr8)*(scale_pyr), "b-2", YRANGE=[0,max(sqrt(pyr8)*scale_pyr)], XRANGE=[0,size(pyr8,/dim)], TITLE='Pyramid WFS', YTITLE='noise propagation coefficient [nm]', XTITLE='mode number', dim=[1000,600],NAME='8 modes')
;p2 = plot(sqrt(pyr19)*(scale_pyr), "r-2", XRANGE=[0,size(pyr19,/dim)],YRANGE=[0,max(sqrt(pyr8)*scale_pyr)], /OVERPLOT, NAME='19 modes')
;leg = legend(TARGET=[p,p2],POSITION=[0.9,0.8])
;;p = plot(W, "r", XRANGE=[0,size(W,/dim)], YTITLE='noise propagation coefficient', XTITLE='index', dim=[1000,600])

print, 'done'



end

