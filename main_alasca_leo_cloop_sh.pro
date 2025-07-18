; run on GPU
detect_dlms, /parallel
set_white_plot

tic

; read parameters file
dir = ROUTINE_DIR() 
savedir = '/home2/bstadler/ALASCA_LEO/'
params_file = 'params_alasca_leo_sh_test.pro'
params = read_params_file(dir+params_file)

; ********************************************
;       -----------------
; ----> KEYWORDS
;       -----------------
; ********************************************
objectDir = 0.
pwfs = 0B
ttCorr4Psf = 0B

do_opt_gain = 0B
do_forg_fact = 0B
do_aberration = 1B

; ********************************************
;       -----------------
; ----> PARAMETERS UPDATE
;       -----------------
; ********************************************
rad2arcsec = 3600.*360./2/!pi
diam = params.main.pixel_pupil*params.main.pixel_pitch

if do_opt_gain then begin
  params.control.opt_dt = max([1,params.detector.dt*1000.])
  if params.main.total_time lt 2*params.control.opt_dt then params.main.total_time = 5*params.control.opt_dt
  params.atmo.pixel_phasescreens = 2*32768L
endif
if do_forg_fact then begin
  nmodes = params.modalrec.nmodes
  ff_min_value = 0.986-0.06*params.seeing.constant
  ff_exp = 4.
  imode = 5
  ff = fltarr(nmodes)
  zern_degree, findgen(nmodes)+2, rd, af
  for i=0,nmodes-1 do  ff[i] = min([1, 1 - (1-ff_min_value) * ( float(rd[i]-rd[imode])/float(max(rd)-rd[imode]) )^ff_exp ])
  ff[0:imode-1] = 1.
  params.control.ff = ff
endif

sn_tag = params.main.instrument_name+'_ps120p0.008_shs10x10_wl589_fv'+strtrim(string(params.sh.sensor_fov,format='(f9.1)'),2)+$
  '_np24_th0.25_ce'
params.slopec.sn_tag = sn_tag
ngs_sn_tag = params.main.instrument_name+'_ps120p0.008_shs10x10_wl1064_fv'+strtrim(string(params.sh_ngs.sensor_fov,format='(f9.1)'),2)+$
  '_np24_th0.25_ce'
params.slopec_IR_ngs.sn_tag = ngs_sn_tag

; physical propagation
params.atmo_down_589 = duplicate_params(params.atmo)
params.atmo_down_589.doFresnel      = 1B
params.atmo_down_589.wavelengthInNm = 589
params.atmo_up_1550 = duplicate_params(params.atmo)
params.atmo_up_1550.doFresnel      = 1B
params.atmo_up_1550.wavelengthInNm = 1550
params.atmo_down_1064 = duplicate_params(params.atmo)
params.atmo_down_1064.doFresnel      = 1B
params.atmo_down_1064.wavelengthInNm = 1064

; DM
params.dm_up = duplicate_params(params.dm)
params.dm_IR_ngs = duplicate_params(params.dm)

; automatic tag update

params = update_alasca_leo_params(params, objectDir=objectDir, pwfs=pwfs)


; ********************************************
;       -------
; ----> OBJECTS
;       -------
; ********************************************
; Initialize housekeeping objects
factory = obj_new('factory',params.main, /GPU)
loop    = factory.get_loop_control()
store   = factory.get_datastore()

; ------------------------- Sources (LGS, SAT) -------------------------
source_down_589_list = list()
source_down_589_list.add, factory.get_source(params.wfs_source)
source_up_1550_list = list()
source_up_1550_list.add, factory.get_source(params.wfs_sat_source)
source_down_1064_list = list()
source_down_1064_list.add, factory.get_source(params.wfs_ngs_source)


; ------------------------- Pupilstop -------------------------
pupstop_up   = factory.get_pupilstop( params.pupil_stop_up )
pupstop_down = factory.get_pupilstop( params.pupil_stop_down )


; ------------------------- Atmosphere and propagation -------------------------
atmo_down_589 = factory.get_atmo_container(source_down_589_list, params.atmo_down_589, $
  params.seeing, params.wind_speed, params.wind_direction) ;down-ward atmosphere (same parameters as up one's)
prop_down_589 = factory.get_atmo_propagation(params.atmo_down_589, source_down_589_list) ;down-ward propagation
  
atmo_up_1550   = factory.get_atmo_container(source_up_1550_list, params.atmo_up_1550, $
  params.seeing, params.wind_speed, params.wind_direction) ;up-ward atmosphere
prop_up_1550   = factory.get_atmo_propagation(params.atmo_up_1550, source_up_1550_list) ;up-ward propagation
  
atmo_down_1064 = factory.get_atmo_container(source_down_1064_list, params.atmo_down_1064, $
    params.seeing, params.wind_speed, params.wind_direction) ;down-ward atmosphere (same parameters as up one's) 
prop_down_1064 = factory.get_atmo_propagation(params.atmo_down_1064, source_down_1064_list) ;down-ward propagation


; ------------------------- Wavefront Sensor -------------------------
; for LGS
sh_lgs = factory.get_sh(params.sh)
ccd_lgs = factory.get_ccd(params.detector)
sc_lgs = factory.get_sh_slopec(params.slopec)
ef_prod_sat = factory.get_ef_product()
; extended object
sh_lgs.kernel_application = 'FFT'
kernelobj = factory.get_kernel_full(sh_lgs, params.launcher, params.seeing, params.zlayer, params.zprofile)

; for space object
sh_ngs = factory.get_sh(params.sh_ngs)
ccd_ngs = factory.get_ccd(params.detector_IR_ngs)
sc_ngs = factory.get_sh_slopec(params.slopec_IR_ngs)
ef_prod_ngs = factory.get_ef_product()


; ------------------------- Reconustructor and control -------------------------
rec       = factory.get_modalrec(params.modalrec)
intc      = factory.get_control(params.control)

rec_ngs     = factory.get_modalrec(params.modalrec_IR_ngs)
intc_ngs    = factory.get_control(params.control_IR_ngs)


; ------------------------- Deformable mirrors -------------------------
dm        = factory.get_dm(params.dm)
dm_ngs   = factory.get_dm(params.dm_IR_ngs)
comm_up   = obj_new('base_value')
comm_down   = obj_new('base_value')
dm_up     = factory.get_dm(params.dm_up)


; ------------------------- PSFs -------------------------
psf_down  = factory.get_psf(params.camera)
psf_sat   = factory.get_psf(params.camera_IR_sat) ; IR NGS/SAT PSF
psf_ngs   = factory.get_psf(params.camera_IR_ngs) ; IR NGS/SAT PSF

sr_pyr    = obj_new('base_value')
sr_pyr_max= obj_new('base_value')

; ------------------------- Modal analysis -------------------------
if params.hasKey('modalanalysis') then begin
  modan_down_589 = factory.get_modalanalysis(params.modalanalysis)
  modan_down_589.wavelengthInNm = 589
  modan_ngs = factory.get_modalanalysis(params.modalanalysis)
  modan_ngs.wavelengthInNm = 1064
  modan_sat = factory.get_modalanalysis(params.modalanalysis)
  modan_sat.wavelengthInNm = 1550
endif

;----NM add spherical aberration ----
if do_aberration then aber = factory.get_disturbance(params.aberration)

; ********************************************
;       -------------------
; ----> OBJECTS CONNECTIONS
;       -------------------
; ********************************************
; ------------------------- Atmosphere and propagation -------------------------
; down-ward
;prop_down_589.mask = mask
atmo_down_589_layers = atmo_down_589.layer_list
for i=0,atmo_down_589_layers.count()-1 do prop_down_589.add_layer_to_layer_list, atmo_down_589_layers[atmo_down_589_layers.count()-1-i]
prop_down_589.add_layer_to_layer_list, dm.out_layer
prop_down_589.add_layer_to_layer_list, pupstop_down

; up-ward
prop_up_1550.isUpwards = 1B
;prop_up_1550.mask = pupstop_down.A
prop_up_1550.add_layer_to_layer_list, pupstop_up
;prop_up_1550.add_layer_to_layer_list, dm_lgs_tt.out_layer
prop_up_1550.add_layer_to_layer_list, dm_up.out_layer
if do_aberration then prop_up_1550.add_layer_to_layer_list,  aber.out_layer
atmo_up_1550_layers   = atmo_up_1550.layer_list
foreach layer,atmo_up_1550_layers do prop_up_1550.add_layer_to_layer_list, layer

; down-ward
atmo_down_1064_layers = atmo_down_1064.layer_list
for i=0,atmo_down_1064_layers.count()-1 do prop_down_1064.add_layer_to_layer_list, atmo_down_1064_layers[atmo_down_1064_layers.count()-1-i]
prop_down_1064.add_layer_to_layer_list, dm.out_layer
prop_down_1064.add_layer_to_layer_list, pupstop_down

; set up timing in atmospheric propagations
atmo_down_589_atmo = atmo_down_589.get('atmo')
atmo_up_1550_atmo = atmo_up_1550.get('atmo')
atmo_down_1064_atmo = atmo_down_1064.get('atmo')
airmass = 1.0/cos(params.main.zenithAngleInDeg/180.*!pi)
time_vect_up = (2*params.wfs_source.height - params.atmo.heights)*airmass / 299792458d
time_vect_down = params.atmo.heights*airmass / 299792458d
print, 'atmo layers time vector up-ward [s]  : ', time_vect_up
print, 'atmo layers time vector down-ward [s]: ', time_vect_down
atmo_down_589_atmo.extra_delta_time = time_vect_down
atmo_up_1550_atmo.extra_delta_time = time_vect_up
atmo_down_1064_atmo.extra_delta_time = time_vect_down


; find position
ao_delay = 2
time_diff = time_vect_up-time_vect_down + ao_delay*params.main.time_step
x_time = time_diff*params.wind_speed.constant*airmass
x_zenith = tan(params.main.zenithAngleInDeg/180.*!pi)*params.atmo.heights-x_time
alpha = atan(x_zenith/params.atmo.heights)*180./!pi
alpha_time = (params.main.zenithAngleInDeg-alpha)*3600.

if ttCorr4Psf then begin
  modan_down_tt = factory.get_modalanalysis( $
    dictionary({type: 'zernike', nmodes: 2, pupil_mask_tag: params.pupil_stop_down.pupil_mask_tag}))
  modan_down_tt.wavelengthInNm = 589
  dm_down_tt = factory.get_dm( $
    dictionary({type: 'zernike', nmodes: 2, pupil_mask_tag: params.pupil_stop_down.pupil_mask_tag, height: 0}))
  ef_prod_down = factory.get_ef_product()
endif


; ------------------------- Wavefront Sensor and Mirrors -------------------------
; for LGS
sh_lgs.kernelobj = kernelobj
lgstt = factory.get_lgstt_container((prop_down_589.pupil_list)[0], params.tt_modalAnalysis, params.tt_dm, res_params=params.lgsttres); idea Roberto - use the TT commands from tt_dm for correcting uplink - TT from LGS
sh_lgs.in_ef = lgstt.out_ef
ccd_lgs.in_i = sh_lgs.out_i
sc_lgs.in_pixels = ccd_lgs.out_pixels
rec.in_slopes = sc_lgs.out_slopes  ; slopes objects
intc.in_delta_comm= rec.out_modes  ; modal residual objects

dm_up.in_command = comm_up
dm.in_command = intc.out_comm  ; modal commands objects

ef_prod_sat.in_ef1 = (prop_up_1550.pupil_list)[0]
ef_prod_sat.in_ef2 = dm_ngs.out_layer ; change this here with commands form tt_dm of LGS??

dm_ngs.in_command = intc_ngs.out_comm
ef_prod_ngs.in_ef1 = (prop_down_1064.pupil_list)[0]
ef_prod_ngs.in_ef2 = dm_ngs.out_layer

; for NGS
sh_ngs.in_ef = ef_prod_ngs.out_ef
ccd_ngs.in_i = sh_ngs.out_i
sc_ngs.in_pixels = ccd_ngs.out_pixels
rec_ngs.in_slopes       = sc_ngs.out_slopes ; slopes objects
intc_ngs.in_delta_comm  = rec_ngs.out_modes ; modal residual objects


; ------------------------- PSFs -------------------------
if ttCorr4Psf then begin
  modan_down_tt.in_ef   = (prop_down_589.pupil_list)[0]
  dm_down_tt.in_command = modan_down_tt.out_modes
  ef_prod_down.in_ef1   = (prop_down_589.pupil_list)[0]
  ef_prod_down.in_ef2   = dm_down_tt.out_layer
  psf_down.in_ef        = ef_prod_down.out_ef
endif else begin
  psf_down.in_ef = (prop_down_589.pupil_list)[0] ; electric field from down-ward propagation
endelse

psf_ngs.in_ef = ef_prod_ngs.out_ef             ; electric field from down-ward propagation (plus TT correction from IR)
psf_sat.in_ef = ef_prod_sat.out_ef             ; electric field from down-ward propagation (plus TT correction from IR)

; ------------------------- Modal analysis -------------------------
if params.hasKey('modalanalysis') then begin
  modan_down_589.in_ef = (prop_down_589.pupil_list)[0]
  modan_ngs.in_ef = ef_prod_ngs.out_ef
  modan_sat.in_ef = ef_prod_sat.out_ef
endif

coeff0 = -1/4.
coeff1 = -1/4.
z2 = zern2phi(params.main.pixel_pupil, 2, mask=mask, no_round_mask=no_round_mask, $
  xsign=xsign, ysign=ysign, rot_angle=rot_angle, verbose=verbose)
phase = coeff0*reform(z2[0,*,*])+coeff1*reform(z2[1,*,*])

psf_sat_ref = calc_psf(phase, congrid(pupstop_up.A,params.main.pixel_pupil,params.main.pixel_pupil), imwidth=spsf, normalize=0B, nocenter=1B, GPU=0B)
psf_ngs_ref = calc_psf(phase, congrid(pupstop_down.A,params.main.pixel_pupil,params.main.pixel_pupil), imwidth=spsf, normalize=0B, nocenter=0B, GPU=0B)

; RADIAL PROFILES for FWHM computation
spsf = params.camera.nd*params.main.pixel_pupil

; ********************************************
;       -----------------
; ----> DATA TO BE SAVED
;       -----------------
; ********************************************
; set store data
;store.add, sc_lgs.out_slopes, name='slopes'
;store.add, rec.out_modes, name='deltaComm'
;store.add, intc.out_comm, name='comm'
;store.add, comm_up, name='comm_up'
;store.add, sc_ngs.out_slopes, name='slopes_ngs'
;store.add, rec_ngs.out_modes, name='deltaComm_ngs'
;store.add, intc_ngs.out_comm, name='comm_ngs'
;store.add, psf_down.out_sr, name='sr_down'
store.add, psf_ngs.out_sr, name='sr_ngs'
store.add, psf_sat.out_sr, name='sr_sat'
;store.add, (prop_up_1550.pupil_list)[0], name='res_up_sat_ef'
;store.add, (prop_down_589.pupil_list)[0], name='res_down_ef'
;store.add, (prop_down_1064.pupil_list)[0], name='res_down_ngs_ef'
;store.add, psf_down.out_psf, name='psf_down'
;store.add, psf_ngs.out_psf, name='psf_ngs'
;store.add, psf_sat.out_psf, name='psf_sat'
store.add, psf_ngs.psf_center, name='psf_ngs_center'
store.add, psf_sat.psf_center, name='psf_sat_center'
;if do_opt_gain then store.add, intc.optgain, name='optgain'
;if params.hasKey('modalanalysis') then begin
;  store.add, modan_down_589.out_modes, name='resMod_down'
;  store.add, modan_ngs.out_modes, name='resMod_ngs'
;  store.add, modan_sat.out_modes, name='resMod_sat'
;endif
; -------------------------

ph_up_disp   = factory.get_phase_display((prop_up_1550.pupil_list)[0])   ; residual phase display for up-ward propagation
ph_up_disp.window = 11
ph_up_disp.title = 'PHASE UP'
ph_up_disp.disp_factor = 2

ph_down_disp = factory.get_phase_display((prop_down_1064.pupil_list)[0]) ; residual phase display for down-ward propagation
ph_down_disp.window = 14
ph_down_disp.title = 'PHASE DOWN'
ph_down_disp.disp_factor = 2



; ********************************************
;       ----
; ----> LOOP
;       ----
; ********************************************
; Build loop
loop.add, atmo_down_589
loop.add, atmo_up_1550
loop.add, atmo_down_1064
loop.add, prop_down_589
if do_aberration then loop.add, aber
loop.add, prop_up_1550
loop.add, prop_down_1064
if obj_valid(lgstt) then loop.add, lgstt
loop.add, sh_lgs
loop.add, ccd_lgs
loop.add, sc_lgs
loop.add, rec
loop.add, intc
loop.add, cheat('temp = c_down.value & c.value=temp & c.generation_time = t', c=comm_up, c_down=intc.out_comm)
loop.add, dm
loop.add, dm_up
if ttCorr4Psf then begin
  loop.add, modan_down_tt
  loop.add, dm_down_tt
  loop.add, ef_prod_down
endif
; -------
; NGS
loop.add, ef_prod_sat
loop.add, ef_prod_ngs
loop.add, sh_ngs
loop.add, ccd_ngs
loop.add, sc_ngs
loop.add, rec_ngs
loop.add, intc_ngs
loop.add, dm_ngs
;loop.add, dm_lgs_tt
; -------
loop.add, psf_down
loop.add, psf_ngs
loop.add, psf_sat
wave_string = strtrim(round(params.camera.wavelengthInNm),2)
wave_string_sat = strtrim(round(params.camera_IR_sat.wavelengthInNm),2)
wave_string_ngs = strtrim(round(params.camera_IR_ngs.wavelengthInNm),2)
loop.add, cheat('print, "SR(@'+wave_string+'nm, down-ward) = "+strtrim(sr1.value,2)', sr1=psf_down.out_sr)
loop.add, cheat('print, "SR(@'+wave_string_sat+'nm, up-ward, SAT) = "+strtrim(sr.value,2)', sr=psf_sat.out_sr)
loop.add, cheat('print, "SR(@'+wave_string_ngs+'nm, down-ward, NGS) = "+strtrim(sr.value,2)', sr=psf_ngs.out_sr)
if params.hasKey('modalanalysis') then begin
  loop.add, modan_down_589
  loop.add, modan_ngs
  loop.add, modan_sat
  fitErr = calc_fitting_error(params.seeing.constant,0.1,500e-9,diam=diam,zenith=params.main.zenithAngleInDeg,/returnErrorInNm)
endif
loop.add, store
;loop.add, ph_up_disp
;loop.add, ph_down_disp


; ********************************************
;       ---
; ----> RUN
;       ---
; ********************************************
; Run simulation loop
loop.run, run_time=params.main.total_time, dt=params.main.time_step, /speed_report


; ********************************************
;       ----------------------------------------
; ----> ADDITIONAL DATA, FINAL PRINTS AND SAVING
;       ----------------------------------------
; ********************************************
;add integrated PSF to store
store.add_array, psf_down.out_int_psf.value, name='int_psf_down'
store.add_array, psf_ngs.out_int_psf.value, name='int_psf_ngs'
store.add_array, psf_sat.out_int_psf.value, name='int_psf_sat'

; saving
tn = store.save_tracknum(dir=savedir, params=params, /nodlm, /nooldformat, /compress, /saveFloat)

; prints
print, 'Mean SE Strehl Ratio (@'+strtrim(params.camera.wavelengthInNm,2)+'nm, down-ward) :', $
  store.mean('sr_down',init=min([50,0.1*params.main.total_time/params.main.time_step]))*100., '%'
print, 'Mean SE Strehl Ratio (@'+strtrim(params.camera_IR_ngs.wavelengthInNm,2)+'nm, down-ward) :', $
  store.mean('sr_ngs',init=min([50,0.1*params.main.total_time/params.main.time_step]))*100., '%'
print, 'Mean SE Strehl Ratio (@'+strtrim(params.camera_IR_sat.wavelengthInNm,2)+'nm, SAT) :', $
  store.mean('sr_sat',init=min([50,0.1*params.main.total_time/params.main.time_step]))*100., '%'

print, 'LE Strehl Ratio (@'+strtrim(params.camera.wavelengthInNm,2)+'nm, down-ward) :', $
  (psf_down.out_int_psf.value/total(psf_down.out_int_psf.value))[spsf/2,spsf/2]/psf_down.out_ref_sr*100., '%'
print, 'LE Strehl Ratio (@'+strtrim(params.camera_IR_ngs.wavelengthInNm,2)+'nm, down-ward) :', $
  (psf_ngs.out_int_psf.value/total(psf_ngs.out_int_psf.value))[spsf/2,spsf/2]/psf_ngs.out_ref_sr*100., '%'
print, 'LE Strehl Ratio (@'+strtrim(params.camera_IR_sat.wavelengthInNm,2)+'nm, SAT) :', $
  (psf_sat.out_int_psf.value/total(psf_sat.out_int_psf.value))[spsf/2,spsf/2]/psf_sat.out_ref_sr*100., '%'
print, 'Scintillation, down-ward :', variance((prop_down_1064.power).ToArray()/mean( (prop_down_1064.power).ToArray()))
  
print, psf_ngs_ref[0,0]
print, psf_sat_ref[0,0]

toc
end
