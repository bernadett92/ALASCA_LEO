detect_dlms, /par

objectDir = 90.
pwfs = 0B
over = 0B

dir = ROUTINE_DIR()
params = (read_params_file(dir+'params_alasca_leo_sh.pro',/expand))[0]
subap_params = read_params_file(dir+'params_alasca_leo_subaps_sh.pro')
sn_params = read_params_file(dir+'params_alasca_leo_sn_sh.pro')
intmat_params = read_params_file(dir+'params_alasca_leo_intmat_sh.pro')

params = update_alasca_leo_params(params, objectDir=objectDir, pwfs=pwfs)

params.remove, 'launcher'
params.remove, 'zlayer'
params.remove, 'zprofile'
params.remove, 'lgsttres'

; ------------------------- LGS TT reconstruction matrix -------------------------
params_tt = duplicate_params(params)
params_tt.wfs_source = params_tt.wfs_source
params_tt.sh = params_tt.sh
params_tt.detector = params_tt.detector
params_tt.slopec = params_tt.slopec
params_tt.modalrec = params_tt.modalrec_tt

intmat_params_tt = duplicate_params(intmat_params)
subap_params_tt = duplicate_params(subap_params)
sn_params_tt = duplicate_params(sn_params)

calib_tt = obj_new('scao_calib', params_tt, GPU=1B, display=1B, disp_factor=disp_factor)
calib_tt.buildSubaps, subapsParams=subap_params_tt, subapsExist=subapsExist_tt
if subapsExist_tt eq 0 or over then calib_tt.runSubaps
calib_tt.buildSlopeNull, slopeNullParams=sn_params_tt, slopeNullExist=slopeNullExist_tt
if slopeNullExist_tt eq 0 or over then calib_tt.runSlopeNull
calib_tt.buildIntmat, intmatParams=intmat_params_tt, intmatExist=intmatExist_tt, recmatExist=recmatExist_tt
if intmatExist_tt eq 0 or over then calib_tt.runIntmat
if recmatExist_tt eq 0 or over then calib_tt.computeRecmat

end
