
detect_dlms, /parallel
seeing = 1.0
seing4calib = 1.0
objectDir = 0.
pwfs = 1B

dir = ROUTINE_DIR()
if seeing eq 2.5 then params_file='params_alasca_leo_daytime.pro' else params_file = 'params_alasca_leo_pyr.pro'
if seeing eq 5.0 then params_file='params_alasca_leo_daytime_5.0.pro'
params = (read_params_file(dir+params_file,/expand))[0]
pupils_params = read_params_file(dir+'params_alasca_leo_pupils.pro')
sn_params = read_params_file(dir+'params_alasca_leo_sn.pro')
intmat_params = read_params_file(dir+'params_alasca_leo_intmat.pro')

over = 0B
params = update_alasca_leo_params(params, objectDir=objectDir, pwfs=pwfs)
params.remove, 'extended_object'

disp_factor=2
calib = obj_new('scao_calib', params, GPU=1B, display=0B, disp_factor=disp_factor)

tic
calib.buildPupils, pupilsParams=pupils_params, pupilsExist=pupilsExist
if pupilsExist eq 0 or over then calib.runPupils
calib.buildSlopeNull, slopeNullParams=sn_params, slopeNullExist=slopeNullExist
if slopeNullExist eq 0 or over then calib.runSlopeNull
toc

end

