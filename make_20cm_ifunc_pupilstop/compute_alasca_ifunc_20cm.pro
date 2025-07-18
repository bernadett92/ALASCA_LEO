IF_DIR = '/home/bstadler/Simulator/PASSATA/mainPASSATA/ALASCA_LEO/make_20cm_ifunc_pupilstop/ALPAO_IF/DM_Calib_2021.12.09_wtt/act_wtt_'
OUTPUT_DIR = '/home/bstadler/passata/ALASCA_LEO/ifunc/'
N_OUTPUT = 24L ; = 120*0.2cm
IFMASK = make_mask(N_OUTPUT,obs=0.0)
BAD_ACTS = [31, 43, 75]

ifunc=read_alpao_if(IF_DIR, ifmask=ifmaskTemp, VERBOSE=VERBOSE, DISPLAY=DISPLAY, MASK_ACTS=MASK_ACTS)

original_nacts = (size(ifunc,/dim))[2]

; Bad acts removal
if keyword_set(BAD_ACTS) then begin
   idx = lonarr(original_nacts)+1
   idx[BAD_ACTS]=0
   use_acts = where(idx)
   ifunc = ifunc[*,*,use_acts]
   nacts = n_elements(use_acts)
endif

ifunc = ifunc[24:270,24:270,*]

; resize
if n_elements(N_OUTPUT) gt 0 then begin
    ifunc_small = make_array(N_OUTPUT,N_OUTPUT,(size(ifunc,/dim))[2],type=size(ifunc,/type))
    mask_big = float(reform(total(abs(ifunc),3)) gt 1e-6)
    idxBig = where(mask_big)
    mask_small = make_mask(N_OUTPUT)
    idxSmall = where(mask_small,nMaskFull)
    for i=0,(size(ifunc,/dim))[2]-1 do begin
        if ~ptr_valid(extrapol_mat) then begin
            extrapol_mat_X = ptr_new(extrapolate_edge_pixel_mat_define(mask_big,/doExt2Pix,/onlyX))
            extrapol_mat_Y = ptr_new(extrapolate_edge_pixel_mat_define(mask_big,/doExt2Pix,/onlyY))
        endif
               
        tempAX = reform(ifunc[*,*,i])
        tempAX = extrapolate_edge_pixel(tempAX, *(extrapol_mat_X))
        tempAY = tempAX
        tempAY = extrapolate_edge_pixel(tempAY, *(extrapol_mat_Y)) 
        temp = tempAY
       
        temp = toccd(temp,N_OUTPUT,GPU=0B)
        ifunc_small[*,*,i] = temp * float(mask_small)
    endfor
    ifunc = ifunc_small
endif

sIfunc = size(ifunc,/dim)
piston_orig = reform(total(ifunc,3))
xx_coord = fltarr(sIfunc[2])
yy_coord = fltarr(sIfunc[2])
p2v = fltarr(sIfunc[2])

; actuators coordinates
make_xy, sIfunc[1], sIfunc[1]/2., xx, yy 
for i=0,sIfunc[2]-1 do begin
    dummy = max(smooth(ifunc_small[*,*,i],2),idx_max)
    xx_coord[i] = xx[idx_max]
    yy_coord[i] = yy[idx_max]
    p2v[i] = max(ifunc[*,*,i]) - min(ifunc[*,*,i])
    ; ------
    maskDistance = 0B
    if maskDistance then begin
        distance_i = sqrt( (xx-xx_coord[i])^2. + (yy-yy_coord[i])^2. )
        sort_distance_i = sort(distance_i)
        temp = ifunc[*,*,i] - mean(ifunc,dim=3)
        nTemp = n_elements(sort_distance_i)-n_elements(sort_distance_i)/8
        temp[sort_distance_i[n_elements(sort_distance_i)/8:*]] *= 1-findgen(nTemp)/nTemp
        ifunc[*,*,i] = temp
    endif
endfor
distance_from_0 = sqrt(xx_coord^2. + yy_coord^2.)

ifunc2d_full = fltarr(nacts,nMaskFull)
for i=0,nacts-1 do ifunc2d_full[i,*] = (ifunc[*,*,i])[idxSmall]

; slaving
idx_master = where(p2v gt 0.82*max(p2v) or distance_from_0 lt 0.95*sIfunc[0]/2., nmaster, complement=idx_slave, ncomplement=nslave)
if nslave gt 0 then begin
    print, nslave, ' actuators are slaves'
    print, nmaster, ' actuators are masters'
    slaving_mat = fltarr(sIfunc[2],nmaster)
    for i=0,nmaster-1 do slaving_mat[idx_master[i],i] = 1.
    ifunc_temp = ifunc[*,*,idx_master]
    for i=0,nslave-1 do begin
        distance = sqrt( (xx_coord[idx_master] - xx_coord[idx_slave[i]])^2. + (yy_coord[idx_master] - yy_coord[idx_slave[i]])^2. ) 
        dSort = sort(distance)
        ifunc_temp[*,*,dSort[0]] += ifunc[*,*,idx_slave[i]]
        slaving_mat[idx_slave[i],dSort[0]] += 1.
    endfor
    ifunc = ifunc_temp    
endif

piston = reform(total(ifunc,3))
sIfunc = size(ifunc,/dim)

set_white_plot
loadct, 3
window, 0, xs=1500,ys=500
image_show, /as, /sh, [piston_orig/max(piston_orig),piston/max(piston),piston_orig/max(piston_orig)-piston/max(piston)]

id_mode_starting = 0L
nRawCol = 10
if nRawCol ge 8 then begin
    xs = 1600
    ys = 1200
endif else begin
    xs = 1200
    ys = 900
endelse

shapeBig = fltarr(nRawCol*sIfunc[0],nRawCol*sIfunc[0])
for i=0,nRawCol-1 do begin
    for j=0,nRawCol-1 do begin
        id_mode = id_mode_starting+i+(nRawCol-1-j)*nRawCol
        if id_mode le nmaster-1 then begin
            shapeBig[i*sIfunc[0]:(i+1)*sIfunc[0]-1,j*sIfunc[0]:(j+1)*sIfunc[0]-1] = ifunc_temp[*,*,id_mode]
        endif
    endfor
endfor
    
window, 1, xs=xs, ys=ys
image_show, /as, /sh, shapeBig, tit='!17Zonal -  mode '+strtrim(round(id_mode_starting),2)+'-'+strtrim(round(id_mode_starting+nRawCol^2.-1),2)
idx_IFMASK = where(IFMASK,nMask)
ifunc2d = fltarr(nmaster,nMask)
for i=0,nmaster-1 do ifunc2d[i,*] = (ifunc[*,*,i])[idx_IFMASK]

; Some default parameters...
diam = 1.0
r0 = 0.5
L0 = 22.0
ZernModes = 2
IFmaxConditionNumber = 1e4
oversampling = 2

klbasis = make_modal_base_from_ifs_fft(IFMASK, diam, ifunc2d, $
                    r0, L0, ZernModes=ZernModes, oversampling=oversampling, M2C=M2C, $
                    IFmaxConditionNumber=IFmaxConditionNumber, VERBOSE=1B, DOUBLE=0B)
IFMASK_big = fltarr(120,120)
IFMASK_big[60-N_OUTPUT/2:60+N_OUTPUT/2-1,60-N_OUTPUT/2:60+N_OUTPUT/2-1] = IFMASK
IFMASK = shift(IFMASK_big,40)
idx_IFMASK = where(IFMASK_big,nMask)
N_OUTPUT = 120
sIfunc = [120,120,sIfunc[2]]

removeMode2 = 1B
if removeMode2 then begin
    idx_good = [0,1,findgen((size(klbasis,/dim))[0]-3)+3]
    klbasis = klbasis[idx_good,*]               
    M2C = M2C[idx_good,*]
    print
    print, "WARNING: mode 2 was removed!"
    print
endif

cov = matmat_multiply(klbasis,do_transpose(klbasis))/nMask
covExtraDiag = cov - diag_matrix(diag_matrix(cov))

id_mode_starting = 0L
nRawCol = 8
if nRawCol ge 8 then begin
    xs = 1600
    ys = 1200
endif else begin
    xs = 1200
    ys = 900
endelse

shapeBig = fltarr(nRawCol*sIfunc[0],nRawCol*sIfunc[0])
for i=0,nRawCol-1 do begin
    for j=0,nRawCol-1 do begin
        id_mode = id_mode_starting+i+(nRawCol-1-j)*nRawCol
        if id_mode le (size(klbasis,/dim))[0]-1 then begin
            shape = fltarr(sIfunc[0],sIfunc[0])
            shape[idx_IFMASK] = klbasis[id_mode,*]
            shapeBig[i*sIfunc[0]:(i+1)*sIfunc[0]-1,j*sIfunc[0]:(j+1)*sIfunc[0]-1] = shape
        endif
    endfor
endfor

loadct, 3
window, 2, xs=800, ys=600
plot_oi, findgen(n_elements(diag_matrix(cov)))+1, diag_matrix(cov), /ynozero, /xst, tit='!17', ytit='covriance matrix diagonal values', xtit='mode number'
window, 3, xs=800, ys=600
image_show, /as, /sh, abs(covExtraDiag)>1e-9, /log, tit='!17abs. of Covariance matrix without diagonal'
window, 4, xs=xs, ys=ys
image_show, /as, /sh, shapeBig, tit='!17fft KL fitting - mode '+strtrim(round(id_mode_starting),2)+'-'+strtrim(round(id_mode_starting+nRawCol^2.-1),2), $
    min=-3, max=3

m2c = matrix_multiply(m2c,slaving_mat,/btr)

; Adjust M2C for the removed bad actuators
nmodes = n_elements(m2c[*,0])
nacts  = n_elements(m2c[0,*])

m2c_full = fltarr(nmodes, original_nacts)
m2c_full[*,use_acts] = m2c
m2c = m2c_full

ifunc2d_full_full = fltarr(original_nacts, (size(ifunc2d_full,/dim))[1])
ifunc2d_full_full[use_acts,*] = ifunc2d_full
ifunc2d_full = ifunc2d_full_full

stop

; *********************
; build m2c object and save it
m2c_obj = obj_new('m2c')
m2c_obj.m2c = M2C
m2c_obj.save, OUTPUT_DIR+'m2c/ALASCA_LEO_20cm_dm_m2c_'+strtrim(nmaster,2)+'masters_'+strtrim((size(klbasis,/dim))[0],2)+'modes.fits'
; *********************
; build modal ifunc object and save it
ifunc = obj_new('ifunc')
ifunc.influence_function = ifunc2d_full
ifunc.mask_inf_func = IFMASK
sxaddpar, hdr, 'INF_FUNC', 'zonalIFs'
sxaddpar, hdr, 'N__MODES', strtrim((size(ifunc2d_full,/dim))[0],2)
sxaddpar, hdr, 'X_Y_SIDE', strtrim(sIfunc[0],2)
ifunc.save, OUTPUT_DIR+'ifunc/ALASCA_LEO_20cm_dm_zonal_ifunc_'+strtrim(sIfunc[0],2)+'pix.fits'
; *********************
; build modal ifunc object and save it
ifunc = obj_new('ifunc')
ifunc.influence_function = klbasis
ifunc.mask_inf_func = ifmask
sxaddpar, hdr, 'INF_FUNC', 'KLsTTF'
sxaddpar, hdr, 'N__MODES', strtrim((size(klbasis,/dim))[0],2)
sxaddpar, hdr, 'X_Y_SIDE', strtrim(sIfunc[0],2)
ifunc.save, OUTPUT_DIR+'ifunc/ALASCA_LEO_20cm_dm_ifunc_'+strtrim(sIfunc[0],2)+'pix_'+strtrim(nmaster,2)+'masters_'+strtrim((size(klbasis,/dim))[0],2)+'modes.fits'
; *********************
; inverse modal ifunc object and save it
ifunc_inv = obj_new('ifunc')
ifunc_inv.influence_function = pseudo_invert(klbasis)
ifunc_inv.mask_inf_func = ifmask
sxaddpar, hdr, 'INF_FUNC', 'KLsTTF'
sxaddpar, hdr, 'N__MODES', strtrim(nmodes,2)
sxaddpar, hdr, 'X_Y_SIDE', strtrim(sIfunc[0],2)
ifunc_inv.save, OUTPUT_DIR+'ifunc/ALASCA_LEO_20cm_dm_ifunc_'+strtrim(sIfunc[0],2)+'pix_'+strtrim(nmaster,2)+'masters_'+strtrim((size(klbasis,/dim))[0],2)+'modes_inv.fits'
end

