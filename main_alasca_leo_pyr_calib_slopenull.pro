detect_dlms, /par

seeing4calib = 1.0
objectDir = 0.
pwfs = 1B

dir = ROUTINE_DIR()
params_file = 'params_alasca_leo_pyr.pro'
params = (read_params_file(dir+params_file,/expand))[0]
pupils_params = read_params_file(dir+'params_alasca_leo_pupils.pro')
sn_params = read_params_file(dir+'params_alasca_leo_sn.pro')
params_file_im = 'params_alasca_leo_intmat.pro'
intmat_params = read_params_file(dir+params_file_im)
params.slopec.use_sn = 0B
if params.hasKey('calib_source') then params.remove, 'calib_source'

nmodes = params.modalrec.nmodes

airmass = 1/cos(params.main.zenithAngleInDeg/180.*!pi)
seeing4calib *= airmass^(3/5.)
seeing4calib = round(seeing4calib*10.)/10.

nframes = 500L
maxFoV = params.pyramid.fov
rad2arcsec = 3600.*360./2/!pi
diam = params.main.pixel_pupil*params.main.pixel_pitch
pixelScalePSF = params.pyramid.wavelengthInNm*1e-9/diam*rad2arcsec/params.pyramid.fft_res
npsf = ceil(maxFoV/pixelScalePSF/2.)*2L
aber_dir = params.main.root_dir+'/data/'
aber_name = 'aber_s'+strtrim(string(seeing4calib,format='(f9.2)'),2)+'asec_fitonly_nm'+strtrim(round(nmodes),2)+'_steps'+strtrim(round(nframes),2)
aber_PSF_name = 'PSF_s'+strtrim(string(seeing4calib,format='(f9.2)'),2)+'asec_fitonly_nm'+strtrim(round(nmodes),2)+'_steps'+$
    strtrim(round(nframes),2)+'_l'+strtrim(round(params.pyramid.wavelengthInNm),2)+'nm_n'+strtrim(round(npsf),2)+'pix'
aber_PSF_bistatic_name = 'PSF_s'+strtrim(string(seeing4calib,format='(f9.2)'),2)+'asec_seeing_w'+strtrim(round(params.pyramid.wavelengthInNm),2)+$
    'nm_n'+strtrim(round(npsf),2)+'pix'

if ~params.detector.hasKey('binning') then binning = 1 else binning = params.detector.binning

sn_tag = params.main.instrument_name+'_ps120p0.008_pyr40x40_wl589_fv'+strtrim(string(params.pyramid.fov,format='(f9.1)'),2)+$
          '_ft3.0_bn'+strtrim(round(binning),2)+'_th0.30a0.30b_s'+$
          strtrim(string(seeing4calib,format='(f9.2)'),2)+'asec_ave_steps'+strtrim(round(nframes),2)
sn_tag += '_extObj'
sn_tag += '_monostatic'
params.slopec.sn_tag = sn_tag

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

params.calib_stat_dist = dictionary({map_tag:aber_name, height:0., dm_npixels:120})

print, 'params.slopec.sn_tag: ', params.slopec.sn_tag

params.detector.readout_noise = 0B
params.detector.photon_noise = 0B
params.detector.background_noise = 0B
params.detector.darkcurrent_noise = 0B
params.detector.excess_noise = 0B

params = update_ALASCA_tags(params,seeing4calib=seeing4calib,bistatic=0B)

pyr_measure_sn, params, params.slopec.sn_tag, calib_params=calib_params, addext=1B, display=0B, GPU=GPU

end