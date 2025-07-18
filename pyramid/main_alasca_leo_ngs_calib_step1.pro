
detect_dlms, /par

objectDir = 0.
pwfs = 1B

dir = ROUTINE_DIR()
params_file = 'params_alasca_leo_pyr.pro'
params = (read_params_file(dir+params_file,/expand))[0]
pupils_params = read_params_file(dir+'params_alasca_leo_pupils.pro')
sn_params = read_params_file(dir+'params_alasca_leo_sn.pro')
intmat_params = read_params_file(dir+'params_alasca_leo_intmat.pro')

pupils_params.calib_source.wavelengthInNm = params.pyramid_IR_ngs.wavelengthInNm
sn_params.calib_source.wavelengthInNm = params.pyramid_IR_ngs.wavelengthInNm
intmat_params.calib_source.wavelengthInNm = params.pyramid_IR_ngs.wavelengthInNm
params.remove, 'extended_object'

over = 0B

params_ngs = update_alasca_leo_params(params, objectDir=objectDir, pwfs=pwfs)

params.wfs_source = params_ngs.wfs_ngs_source
params.pyramid = params_ngs.pyramid_IR_ngs
params.detector = params_ngs.detector_IR_ngs
params.slopec = params_ngs.slopec_IR_ngs
params.modalrec = params_ngs.modalrec_IR_ngs

disp_factor=2
calib = obj_new('scao_calib', params, GPU=1B, display=1B, disp_factor=disp_factor)

tic
calib.buildPupils, pupilsParams=pupils_params, pupilsExist=pupilsExist
if pupilsExist eq 0 or over then calib.runPupils
calib.buildSlopeNull, slopeNullParams=sn_params, slopeNullExist=slopeNullExist
if slopeNullExist eq 0 or over then calib.runSlopeNull
toc

end
