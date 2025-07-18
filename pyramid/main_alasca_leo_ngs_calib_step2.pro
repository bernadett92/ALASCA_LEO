detect_dlms, /par

seeing4calib = 2.5
objectDir = 90.
pwfs = 1B

dir = ROUTINE_DIR()
params_file = 'params_alasca_leo_pyr.pro'
params = (read_params_file(dir+params_file,/expand))[0]
pupils_params = read_params_file(dir+'params_alasca_leo_pupils.pro')
sn_params = read_params_file(dir+'params_alasca_leo_sn.pro')
params_file_im = 'params_alasca_leo_intmat.pro'
intmat_params = read_params_file(dir+params_file_im)

over = 0B
display = 0B

pupils_params.calib_source.wavelengthInNm = params.pyramid_IR_ngs.wavelengthInNm
sn_params.calib_source.wavelengthInNm = params.pyramid_IR_ngs.wavelengthInNm
intmat_params.calib_source.wavelengthInNm = params.pyramid_IR_ngs.wavelengthInNm

params.remove, 'extended_object'

nmodes = params.modalrec_IR_ngs.nmodes
airmass = 1/cos(params.main.zenithAngleInDeg/180.*!pi)
seeing4calib *= airmass^(3/5.)
seeing4calib = round(seeing4calib*10.)/10.
nframes = 100L 
rad2arcsec = 3600.*360./2/!pi
diam = params.main.pixel_pupil*params.main.pixel_pitch
pixelScalePSF = params.pyramid_IR_ngs.wavelengthInNm*1e-9/diam*rad2arcsec/params.pyramid_IR_ngs.fft_res
aber_dir = params.main.root_dir+'/data/'
aber_name = 'aber_s'+strtrim(string(seeing4calib,format='(f9.2)'),2)+'asec_fitonly_nm'+strtrim(round(nmodes),2)+'_steps'+strtrim(round(nframes),2)

if ~file_test(aber_dir+aber_name+'.fits') then begin

    L0 = params.atmo.L0
    directory = params.main.root_dir+'/phasescreens/'
    ifunc_directory = params.main.root_dir+'/ifunc/'
    ifunc_tag = 'CaNaPy_30cm_dm_ifunc_120pix_68masters_66modes'
    seed = 2
    pixel_square_phasescreens = 8192l
    masked = 1B
    imwidth = params.pyramid_IR_ngs.fft_res*params.main.pixel_pupil
    doPsf = 1B

    pupil_mask = !null 

    cube = calc_corr_phase_cube(seeing4calib, params.atmo.L0, params.main.pixel_pitch, params.main.pixel_pupil, nmodes, $
        directory, ifunc_directory, ifunc_tag, nframes, seed=seed, $
        pixel_square_phasescreens=pixel_square_phasescreens, masked=masked, $
        doPsf=doPsf, psf_lambda=params.pyramid_IR_ngs.wavelengthInNm, imwidth=imwidth, $
        outPsf=outPsf, pupil_mask=pupil_mask)
   
    writefits, aber_dir+aber_name+'.fits', cube
    print, aber_dir+aber_name+'.fits SAVED!'
    wait, 1

endif

intmat_params.calib_disturb_generator.ncycles = nframes
intmat_params.calib_stat_dist = dictionary({map_tag: aber_name, map_cycle: 'push-pull', dm_npixels: params.main.pixel_pupil, height : 0})

params_ngs =  update_alasca_leo_params(params, objectDir=objectDir, pwfs=pwfs)

params.wfs_source = params_ngs.wfs_ngs_source
params.pyramid = params_ngs.pyramid_IR_ngs
params.detector = params_ngs.detector_IR_ngs
params.slopec = params_ngs.slopec_IR_ngs
params.modalrec = params_ngs.modalrec_IR_ngs

disp_factor=2
calib = obj_new('scao_calib', params, GPU=1B, display=display, disp_factor=disp_factor)
tic

calib.buildIntmat, intmatParams=intmat_params, intmatExist=intmatExist, recmatExist=recmatExist
if intmatExist eq 0 or over then calib.runIntmat
if recmatExist eq 0 or over then calib.computeRecmat
toc

end
