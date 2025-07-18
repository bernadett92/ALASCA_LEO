; run on GPU
detect_dlms, /parallel
set_white_plot

tic

; read parameters file
dir = ROUTINE_DIR()
savedir = '/home2/bstadler/ALASCA_LEO/'
params_file = 'params_alasca_leo_woWFS.pro'
params = read_params_file(dir+params_file)


; Notes:
; ---------------------------------------------------------------------------------------------------------------------------------
; - Without WFS and DM: 
;         Conjugatd phase directly applied without using a DM, no delay
;         Downlink: Very good results -> no delay, not affected by wind speed
;         Uplink: PAA problem
;  - Without WFS, with DM: 
;         Phase projected on the DM (modal coordinates), integrator control, loop delay, only WFS and reconstruction skipped 
;         Downlink: Worse compared to other case -> delay in combination with wind speed
;         Uplink: PAA problem
useDM = 1B


; ********************************************
;       -----------------
; ----> PARAMETERS UPDATE
;       -----------------
; ********************************************
rad2arcsec = 3600.*360./2/!pi
diam = params.main.pixel_pupil*params.main.pixel_pitch

; physical propagation
params.atmo_up_1550 = duplicate_params(params.atmo)
params.atmo_up_1550.doFresnel      = 1B
params.atmo_up_1550.wavelengthInNm = 1550
params.atmo_down_1064 = duplicate_params(params.atmo)
params.atmo_down_1064.doFresnel      = 1B
params.atmo_down_1064.wavelengthInNm = 1064

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
source_up_1550_list = list()
source_up_1550_list.add, factory.get_source(params.wfs_sat_source)
source_down_1064_list = list()
source_down_1064_list.add, factory.get_source(params.wfs_ngs_source)

; ------------------------- Pupilstop -------------------------
pupstop_up   = factory.get_pupilstop( params.pupil_stop_up )
pupstop_down = factory.get_pupilstop( params.pupil_stop_down )

; ------------------------- Atmosphere and propagation -------------------------
atmo_up_1550   = factory.get_atmo_container(source_up_1550_list, params.atmo_up_1550, $
  params.seeing, params.wind_speed, params.wind_direction) ;up-ward atmosphere
prop_up_1550   = factory.get_atmo_propagation(params.atmo_up_1550, source_up_1550_list) ;up-ward propagation

atmo_down_1064 = factory.get_atmo_container(source_down_1064_list, params.atmo_down_1064, $
  params.seeing, params.wind_speed, params.wind_direction) ;down-ward atmosphere (same parameters as up one's)
  
prop_down_1064 = factory.get_atmo_propagation(params.atmo_down_1064, source_down_1064_list) ;down-ward propagation

prop_down_1064_woDM = factory.get_atmo_propagation(params.atmo_down_1064, source_down_1064_list) ;down-ward propagation


; ------------------------- DM -------------------------
intc      = factory.get_control(params.control)
dm = factory.get_dm(params.dm)
comms   = obj_new('base_value')
modan_ngs = factory.get_modalanalysis(params.modalanalysis)
modan_ngs.wavelengthInNm = 1064

; ------------------------- PSFs -------------------------
psf_sat   = factory.get_psf(params.camera_IR_sat) ; IR NGS/SAT PSF
psf_ngs   = factory.get_psf(params.camera_IR_ngs) ; IR NGS/SAT PSF


; ********************************************d
;       -------------------
; ----> OBJECTS CONNECTIONS
;       -------------------
; ********************************************
down_conjugate = obj_new( 'LAYER', params.main.pixel_pupil, params.main.pixel_pupil, params.main.pixel_pitch, 0, GPU=1B, PRECISION=0B)
down_conjugate.phaseInNm = fltarr(params.main.pixel_pupil,params.main.pixel_pupil)
down_conjugate.A = fltarr(params.main.pixel_pupil,params.main.pixel_pupil)

; up-ward
prop_up_1550.isUpwards = 1B
prop_up_1550.mask = pupstop_down.A
atmo_up_1550_layers   = atmo_up_1550.layer_list
prop_up_1550.add_layer_to_layer_list, pupstop_up
;prop_up_1550.add_layer_to_layer_list, dm_lgs_tt.out_layer
if useDM then prop_up_1550.add_layer_to_layer_list, dm.out_layer $
  else prop_up_1550.add_layer_to_layer_list, down_conjugate
foreach layer,atmo_up_1550_layers do prop_up_1550.add_layer_to_layer_list, layer

; down-ward
atmo_down_1064_layers = atmo_down_1064.layer_list
for i=0,atmo_down_1064_layers.count()-1 do prop_down_1064.add_layer_to_layer_list, atmo_down_1064_layers[atmo_down_1064_layers.count()-1-i]
if useDM then prop_down_1064.add_layer_to_layer_list, dm.out_layer $
  else prop_down_1064.add_layer_to_layer_list, down_conjugate
prop_down_1064.add_layer_to_layer_list, pupstop_down

for i=0,atmo_down_1064_layers.count()-1 do prop_down_1064_woDM.add_layer_to_layer_list, atmo_down_1064_layers[atmo_down_1064_layers.count()-1-i]
prop_down_1064_woDM.add_layer_to_layer_list, pupstop_down

; set up timing in atmospheric propagations
atmo_up_1550_atmo = atmo_up_1550.get('atmo')
atmo_down_1064_atmo = atmo_down_1064.get('atmo')
airmass = 1.0/cos(params.main.zenithAngleInDeg/180.*!pi)
time_vect_up = (2*params.wfs_source.height - params.atmo.heights)*airmass / 299792458d
time_vect_down = params.atmo.heights*airmass / 299792458d
print, 'atmo layers time vector up-ward [s]  : ', time_vect_up
print, 'atmo layers time vector down-ward [s]: ', time_vect_down
atmo_up_1550_atmo.extra_delta_time = time_vect_up
atmo_down_1064_atmo.extra_delta_time = time_vect_down

; find position
time_diff = time_vect_up-time_vect_down
x_time = time_diff*params.wind_speed.constant*airmass
x_zenith = tan(params.main.zenithAngleInDeg/180.*!pi)*params.atmo.heights-x_time
alpha = atan(x_zenith/params.atmo.heights)*180./!pi
alpha_time = (params.main.zenithAngleInDeg-alpha)*3600.

if useDM then modan_ngs.in_ef = (prop_down_1064.pupil_list)[0] $
  else modan_ngs.in_ef = (prop_down_1064_woDM.pupil_list)[0]

intc.in_delta_comm = modan_ngs.out_modes  ; modal residual objects
dm.in_command = intc.out_comm

; ------------------------- PSFs -------------------------
psf_ngs.in_ef = (prop_down_1064.pupil_list)[0]
psf_sat.in_ef = (prop_up_1550.pupil_list)[0]


; ********************************************
;       -----------------
; ----> DATA TO BE SAVED
;       -----------------
; ********************************************
; set store data
store.add, psf_ngs.out_sr, name='sr_ngs'
store.add, psf_sat.out_sr, name='sr_sat'
store.add, psf_ngs.out_psf, name='psf_ngs'
store.add, psf_sat.out_psf, name='psf_sat'
; -------------------------

ph_up_disp   = factory.get_phase_display((prop_up_1550.pupil_list)[0])   ; residual phase display for up-ward propagation
ph_up_disp.window = 2
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
loop.add, atmo_up_1550
loop.add, atmo_down_1064
loop.add, prop_down_1064_woDM
if not useDM then loop.add, cheat('a.phaseInNm = -b.phaseInNm & a.A = b.A & a.generation_time = b.generation_time',a=down_conjugate, b=(prop_down_1064_woDM.pupil_list)[0])
loop.add, prop_down_1064
loop.add, prop_up_1550
loop.add, modan_ngs
loop.add, intc
loop.add, dm
loop.add, psf_ngs
loop.add, psf_sat
wave_string_sat = strtrim(round(params.camera_IR_sat.wavelengthInNm),2)
wave_string_ngs = strtrim(round(params.camera_IR_ngs.wavelengthInNm),2)
loop.add, cheat('print, "SR(@'+wave_string_sat+'nm, up-ward, SAT) = "+strtrim(sr.value,2)', sr=psf_sat.out_sr)
loop.add, cheat('print, "SR(@'+wave_string_ngs+'nm, down-ward, NGS) = "+strtrim(sr.value,2)', sr=psf_ngs.out_sr)
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
store.add_array, psf_ngs.out_int_psf.value, name='int_psf_ngs'
store.add_array, psf_sat.out_int_psf.value, name='int_psf_sat'
store.add_array,  (prop_down_1064.power).ToArray(), name='power'


; saving
tn = store.save_tracknum(dir=savedir, params=params, /nodlm, /nooldformat, /compress, /saveFloat)

; RADIAL PROFILES for FWHM computation
spsf = params.camera_IR_sat.nd*params.main.pixel_pupil
; prints
print, 'Mean SE Strehl Ratio (@'+strtrim(params.camera_IR_ngs.wavelengthInNm,2)+'nm, down-ward) :', $
  store.mean('sr_ngs',init=min([50,0.1*params.main.total_time/params.main.time_step]))*100., '%'
print, 'Mean SE Strehl Ratio (@'+strtrim(params.camera_IR_sat.wavelengthInNm,2)+'nm, SAT) :', $
  store.mean('sr_sat',init=min([50,0.1*params.main.total_time/params.main.time_step]))*100., '%'

print, 'LE Strehl Ratio (@'+strtrim(params.camera_IR_ngs.wavelengthInNm,2)+'nm, down-ward) :', $
  (psf_ngs.out_int_psf.value/total(psf_ngs.out_int_psf.value))[spsf/2,spsf/2]/psf_ngs.out_ref_sr*100., '%'
print, 'LE Strehl Ratio (@'+strtrim(params.camera_IR_sat.wavelengthInNm,2)+'nm, SAT) :', $
  (psf_sat.out_int_psf.value/total(psf_sat.out_int_psf.value))[spsf/2,spsf/2]/psf_sat.out_ref_sr*100., '%'
print, 'Scintillation, down-ward :', variance((prop_down_1064_woDM.power).ToArray()/mean( (prop_down_1064_woDM.power).ToArray()))

toc
end
