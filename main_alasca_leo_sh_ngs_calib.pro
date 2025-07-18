detect_dlms, /par
objectDir = 0.
pwfs = 0B
over = 0B
tic

dir = ROUTINE_DIR()
params_ngs = (read_params_file(dir+'params_alasca_leo_sh_test.pro',/expand))[0]
params_ngs = update_alasca_leo_params(params_ngs, objectDir=objectDir, pwfs=pwfs)

params_ngs.remove, 'launcher'
params_ngs.remove, 'zlayer'
params_ngs.remove, 'zprofile'
params_ngs.remove, 'lgsttres'

params_ngs.wfs_source = params_ngs.wfs_ngs_source
params_ngs.sh = params_ngs.sh_ngs
params_ngs.detector = params_ngs.detector_IR_ngs
params_ngs.slopec = params_ngs.slopec_IR_ngs
params_ngs.modalrec = params_ngs.modalrec_IR_ngs

subap_params_ngs = read_params_file(dir+'params_alasca_leo_subaps_sh.pro')
sn_params_ngs = read_params_file(dir+'params_alasca_leo_sn_sh_ngs.pro')
intmat_params_ngs = read_params_file(dir+'params_alasca_leo_intmat_sh_ngs.pro')

calib_ngs = obj_new('scao_calib', params_ngs, GPU=1B, display=1B, disp_factor=disp_factor)
calib_ngs.buildSubaps, subapsParams=subap_params_ngs, subapsExist=subapsExist
if subapsExist eq 0 or over then calib_ngs.runSubaps
calib_ngs.buildSlopeNull, slopeNullParams=sn_params_ngs, slopeNullExist=slopeNullExist
if slopeNullExist eq 0 or over then calib_ngs.runSlopeNull
calib_ngs.buildIntmat, intmatParams=intmat_params_ngs, intmatExist=intmatExist, recmatExist=recmatExist
if intmatExist eq 0 or over then calib_ngs.runIntmat
if recmatExist eq 0 or over then calib_ngs.computeRecmat
toc

end

