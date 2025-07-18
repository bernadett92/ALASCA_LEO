detect_dlms, /par

objectDir = 0.
pwfs = 0B
over = 0B

dir = ROUTINE_DIR()
params = (read_params_file(dir+'params_alasca_leo_sh_test.pro',/expand))[0]
subap_params = read_params_file(dir+'params_alasca_leo_subaps_sh.pro')
sn_params = read_params_file(dir+'params_alasca_leo_sn_sh.pro')
intmat_params = read_params_file(dir+'params_alasca_leo_intmat_sh.pro')
params = update_alasca_leo_params(params, objectDir=objectDir, pwfs=pwfs)

params.remove, 'launcher'
params.remove, 'zlayer'
params.remove, 'zprofile'
params.remove, 'lgsttres'

; ------------------------- LGS TT reconstruction matrices ------------------------- 
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
calib_tt.buildSubaps, subapsParams=subap_params, subapsExist=subapsExist_tt
if subapsExist_tt eq 0 or over then calib_tt.runSubaps
calib_tt.buildSlopeNull, slopeNullParams=sn_params, slopeNullExist=slopeNullExist_tt
if slopeNullExist_tt eq 0 or over then calib_tt.runSlopeNull
calib_tt.buildIntmat, intmatParams=intmat_params, intmatExist=intmatExist_tt, recmatExist=recmatExist_tt
if intmatExist_tt eq 0 or over then calib_tt.runIntmat
if recmatExist_tt eq 0 or over then calib_tt.computeRecmat

;params_tt_jitter = duplicate_params(params)
;params_tt_jitter.wfs_source = params_tt_jitter.wfs_source
;params_tt_jitter.sh = params_tt_jitter.sh_tt
;params_tt_jitter.detector = params_tt_jitter.detector_tt
;params_tt_jitter.slopec = params_tt_jitter.slopec_tt
;params_tt_jitter.modalrec = params_tt_jitter.modalrec_tt_jitter
;
;intmat_params_tt_jitter = duplicate_params(intmat_params)
;subap_params_tt_jitter = duplicate_params(subap_params)
;sn_params_tt_jitter = duplicate_params(sn_params)
;
;calib_tt_jitter = obj_new('scao_calib', params_tt_jitter, GPU=1B, display=1B, disp_factor=disp_factor)
;calib_tt_jitter.buildSubaps, subapsParams=subap_params_tt_jitter, subapsExist=subapsExist_tt_jiiter
;if subapsExist_tt_jiiter eq 0 or over then calib_tt_jitter.runSubaps
;calib_tt_jitter.buildSlopeNull, slopeNullParams=sn_params_tt_jitter, slopeNullExist=slopeNullExist_tt_jiiter
;if slopeNullExist_tt_jiiter eq 0 or over then calib_tt_jitter.runSlopeNull
;calib_tt_jitter.buildIntmat, intmatParams=intmat_params_tt_jitter, intmatExist=intmatExist_tt_jiiter, recmatExist=recmatExist_tt_jiiter
;if intmatExist_tt_jiiter eq 0 or over then calib_tt_jitter.runIntmat
;if recmatExist_tt_jiiter eq 0 or over then calib_tt_jitter.computeRecmat

; ------------------------- LGS reconstruction matrices -------------------------
params_lgs = duplicate_params(params)
params_lgs.wfs_source = params_lgs.wfs_source
params_lgs.sh = params_lgs.sh
params_lgs.detector = params_lgs.detector
params_lgs.slopec = params_lgs.slopec

intmat_params_lgs = duplicate_params(intmat_params)
subap_params_lgs = duplicate_params(subap_params)
sn_params_lgs = duplicate_params(sn_params)

calib_lgs = obj_new('scao_calib', params_lgs, GPU=1B, display=1B, disp_factor=disp_factor)
calib_lgs.buildSubaps, subapsParams=subap_params_lgs, subapsExist=subapsExist
if subapsExist eq 0 or over then calib_lgs.runSubaps
calib_lgs.buildSlopeNull, slopeNullParams=sn_params_lgs, slopeNullExist=slopeNullExist
if slopeNullExist eq 0 or over then calib_lgs.runSlopeNull
calib_lgs.buildIntmat, intmatParams=intmat_params_lgs, intmatExist=intmatExist, recmatExist=recmatExist
if intmatExist eq 0 or over then calib_lgs.runIntmat
if recmatExist eq 0 or over then calib_lgs.computeRecmat

end
