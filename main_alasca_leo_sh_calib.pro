detect_dlms, /par

objectDir = 0.
pwfs = 0B
over = 0B

dir = ROUTINE_DIR()
params = (read_params_file(dir+'params_alasca_leo_sh_240.pro',/expand))[0]
subap_params = read_params_file(dir+'params_alasca_leo_subaps_sh.pro')
sn_params = read_params_file(dir+'params_alasca_leo_sn_sh.pro')
intmat_params = read_params_file(dir+'params_alasca_leo_intmat_sh.pro')

;params = update_alasca_leo_params(params, objectDir=objectDir, pwfs=pwfs)
params.seeing.constant = 1.0

params.remove, 'launcher'
params.remove, 'zlayer'
params.remove, 'zprofile'
params.remove, 'lgsttres'

; ------------------------- LGS reconstruction matrices ------------------------- 
params_lgs = duplicate_params(params)
params_lgs.wfs_source = params_lgs.wfs_source
params_lgs.sh = params_lgs.sh_lgs
params_lgs.detector = params_lgs.detector
params_lgs.slopec = params_lgs.slopec_lgs
params_lgs.modalrec = params_lgs.modalrec

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
