
detect_dlms, /parallel
seeing4calib = 2.5
objectDir = 90.
pwfs = 1B

dir = ROUTINE_DIR()
params_file='params_alasca_leo_pyr.pro'
params = (read_params_file(dir+params_file,/expand))[0]
pupils_params = read_params_file(dir+'params_alasca_leo_pupils.pro')
sn_params = read_params_file(dir+'params_alasca_leo_sn.pro')
params_file_im = 'params_alasca_leo_intmat.pro'
intmat_params = read_params_file(dir+params_file_im)

over = 0B
display = 0B

nmodes = params.modalrec.nmodes
airmass = 1/cos(params.main.zenithAngleInDeg/180.*!pi)
seeing4calib *= airmass^(3/5.)
seeing4calib = round(seeing4calib*10.)/10.

print, '--> seeing4calib: ', seeing4calib

nframes = 10L
maxFoV = params.pyramid.fov
rad2arcsec = 3600.*360./2/!pi
diam = params.main.pixel_pupil*params.main.pixel_pitch
pixelScalePSF = params.pyramid.wavelengthInNm*1e-9/diam*rad2arcsec/params.pyramid.fft_res
npsf = ceil(maxFoV/pixelScalePSF/2.)*2L
aber_dir = params.main.root_dir+'/data/'
aber_name = 'aber_s'+strtrim(string(seeing4calib,format='(f9.2)'),2)+'asec_fitonly_nm'+strtrim(round(nmodes),2)+'_steps'+strtrim(round(nframes),2)
aber_PSF_name = 'PSF_s'+strtrim(string(seeing4calib,format='(f9.2)'),2)+'asec_fitonly_nm'+strtrim(round(nmodes),2)+'_steps'+$
    strtrim(round(nframes),2)+'_l'+strtrim(round(params.pyramid.wavelengthInNm),2)+'nm_n'+strtrim(round(npsf),2)+'pix'

if ~file_test(aber_dir+aber_name+'.fits') or ~file_test(aber_dir+aber_PSF_name+'.fits') then begin

    L0 = params.atmo.L0
    directory = params.main.root_dir+'/phasescreens/'
    ifunc_directory = params.main.root_dir+'/ifunc/'
    ifunc_tag = 'CaNaPy_30cm_dm_ifunc_120pix_68masters_66modes'
    seed = 2
    pixel_square_phasescreens = 8192l
    masked = 1B
    imwidth = params.pyramid.fft_res*params.main.pixel_pupil
    doPsf = 1B

    pupil_mask = !null 
    
    cube = calc_corr_phase_cube(seeing4calib, params.atmo.L0, params.main.pixel_pitch, params.main.pixel_pupil, nmodes, $
        directory, ifunc_directory, ifunc_tag, nframes, seed=seed, $
        pixel_square_phasescreens=pixel_square_phasescreens, masked=masked, $
        doPsf=doPsf, psf_lambda=params.pyramid.wavelengthInNm, imwidth=imwidth, $
        outPsf=outPsf, pupil_mask=pupil_mask)
   
    sOutPsf = size(outPsf,/dim)
    outPsf = outPsf[sOutPsf[0]/2-npsf/2:sOutPsf[0]/2+npsf/2-1,sOutPsf[1]/2-npsf/2:sOutPsf[1]/2+npsf/2-1]

    writefits, aber_dir+aber_name+'.fits', cube
    writefits, aber_dir+aber_PSF_name+'.fits', outPsf
    print, aber_dir+aber_name+'.fits SAVED!'
    print, aber_dir+aber_PSF_name+'.fits SAVED!'
    wait, 1

endif
params.extended_object.PSF_tag = aber_PSF_name
rad2arcsec = 3600.*360./2/!pi
diam = params.main.pixel_pupil*params.main.pixel_pitch
params.extended_object.pixelScalePSF = params.pyramid.wavelengthInNm*1e-9/diam*rad2arcsec/params.pyramid.fft_res

intmat_params.calib_disturb_generator.ncycles = nframes
intmat_params.calib_stat_dist = dictionary({map_tag: aber_name, map_cycle: 'push-pull', dm_npixels: params.main.pixel_pupil, height : 0})

params = update_alasca_leo_params(params, objectDir=objectDir, pwfs=pwfs)

disp_factor=2
calib = obj_new('scao_calib', params, GPU=1B, display=display, disp_factor=disp_factor)
tic

calib.buildIntmat, intmatParams=intmat_params, intmatExist=intmatExist, recmatExist=recmatExist
if intmatExist eq 0 or over then calib.runIntmat
if recmatExist eq 0 or over then calib.computeRecmat
toc

end
