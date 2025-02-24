{main,
 root_dir:          '/home/bstadler/passata/ALASCA_LEO/calibration_sh/', ; Root directory for calibration manager
 instrument_name:   'ALASCA_LEO'               ; this is 'CaNaPy' to be compatible with the CaNaPy calibration of 589nm path
 pixel_pupil:       120,                   ; Linear dimension of pupil phase array
 pixel_pitch:       0.008333,              ; [m] Pitch of the pupil phase array
 total_time:        1.0000d,                ; [s] Total simulation running time
 time_step:         0.0005d,                ; [s] Simulation time step
 zenithAngleInDeg:  60.0,
 precision:         0,
 verbose:           0B
}

; --------------------------- Deformable mirrors ---------------------------
{DM,
    ifunc_tag  : 'CaNaPy_30cm_dm_zonal_ifunc_120pix'        ; modes-to-commands matrix tag
    m2c_tag    : 'CaNaPy_30cm_dm_m2c_68masters_66modes'     ; DM command-to-phase influence matrix tag
    nmodes     : 19
    height     : 0                           ; DM height [m]
}
{DM_UP,
 type:              'zernike',             ; modes type
 nmodes:            2,                     ; number of modes
 npixels:           120,                   ; linear dimension of DM phase array
 obsratio:          0.0,                   ; obstruction dimension ratio w.r.t. diameter
 height:            0                      ; DM height [m]
}


; --------------------------- Extended object LGS ---------------------------
; spot elongation
{zlayer,
  func_type :    'SIN'
  constant_tag : 'layerHeights'  ; [m] layer heights
  zfocus :       90000.          ; [m] focus value used as reference
  theta :        [0,0]           ; [arcsec] tip/tilt value used as reference
}
{zprofile,
  func_type :    'SIN'
  constant_tag : 'layerIntensities' ; [m] layer intensity, tag name of the fits file stored in params.main.root_dir+/data/
}
{launcher,             ; LGS launcher position
  position : [0.1,0.1,0]   ; position [x,y,z] meters (x eq y position) original position: x=1250mm, y=5350mm, z=96000mm (x/y angle 13.15degree)
  sizeArcsec: 0.0      ; updated inside main file
}


; --------------------------- Sources ---------------------------
; LGS for HO correction
{wfs_source,
  polar_coordinate:  [0.0, 0.0],              ; [arcsec, degrees] source polar coordinates
  magnitude:         6.2,                     ; source magnitude. Uses n_phot.pro
  zeroPoint:         390.d-10                 ; Uses n_phot.pro
  wavelengthInNm:    589,                     ; [nm] wavelength
  height:            90000.                   ; [m] LGS source altitude
}
; Laser ranging uplink
{wfs_sat_source,
  polar_coordinate:  [0.0, 0.0],              ; [arcsec, degrees] source polar coordinates
  magnitude:         0,                       ; source magnitude. Uses n_phot.pro
  wavelengthInNm:    1055,                    ; [nm] wavelength
  height:            400e3                    ; [m] source altitude
}
; Laser ranging downlink
{wfs_ngs_source,
  polar_coordinate:  [0.0, 0.0],              ; [arcsec, degrees] source polar coordinates
  magnitude:         0,                       ; source magnitude. Uses n_phot.pro
  wavelengthInNm:    1064,                    ; [nm] wavelength
  height:            400e3                    ; [m] source altitude
}


; --------------------------- Wavefront sensors ---------------------------
{sh_lgs,
  wavelengthInNm:     589,
  sensor_fov:         10.0,                
  sensor_pxscale:     10.0/24.,            ; [arcsec] WFS pixel FoV (the WFS Fov is sensor_pxscale*sensor_npx)
  sensor_npx:         24,                 ; no. pixels on the side of the subaperture
  subap_on_diameter:  10,                 ; no. pixels on the diameter of the pupil
  energy_th:          0.5,                ; energy threshold for including subapertures (relative)
  subapdata_tag:      'auto'   
}
{slopec_lgs,
  subapdata_tag: 'auto',                  ; tag of the sub-apertures used by the WFS
  sn_tag:        'auto',                  ; tag of the reference slope vector
  thr_value    : 0.0,                     ; [ph] threshold value used in the slope computation (only pixels with photons above this value used in computation)
  weightedPixRad: 1B,
}
{sh_ngs,
  wavelengthInNm:    1064,
  sensor_fov:        10.0,
  sensor_pxscale:    10.0/24.,            ; [arcsec] WFS pixel FoV (the WFS Fov is sensor_pxscale*sensor_npx)
  sensor_npx:        24,                  ; no. pixels on the side of the subaperture
  subap_on_diameter: 10,                   ; no. pixels on the diameter of the pupil
  energy_th:         0.5,                 ; energy threshold for including subapertures (relative)
  subapdata_tag:     'auto'   
}
{slopec_ngs,
  subapdata_tag: 'auto',                  ; tag of the sub-apertures used by the WFS
  sn_tag       : 'auto',                  ; tag of the reference slope vector
  thr_value    :  0.0,                     ; [ph] threshold value used in the slope computation
  weightedPixRad: 0B,
}

; --------------------------- Detectors ---------------------------
{detector,
   size:              [240,240],             ; Detector size in pixels
   dt:                0.0005d,                ; [s] Detector integration time
   bandw:             10,                    ; [nm] Sensor bandwidth
   photon_noise:      1b,                    ; activate photon noise
   readout_noise:     1b,                    ; activate readout noise
   readout_level:     0.8,                   ; readout noise in [e-/pix/frame]
   quantum_eff:       1.0,                    ; quantum efficiency * total transmission
   verbose:           0b
}
{detector_IR_ngs,
  size:              [240,240],             ; Detector size in pixels
  dt:                0.0005d,               ; [s] Detector integration time
  bandw:             300,                   ; [nm] Sensor bandwidth
  photon_noise:      1b,                    ; activate photon noise
  readout_noise:     1b,                    ; activate readout noise
  readout_level:     0.8,                   ; readout noise in [e-/pix/frame]
  quantum_eff:       1.0                    ; quantum efficiency * total transmission
}


; --------------------------- Cameras ---------------------------
{camera,
  wavelengthInNm:    589,                 ; [nm] Imaging wavelength
  nd:                3,                    ; padding coefficient for PSF computation
}
{camera_IR_sat,
  wavelengthInNm:    1055,                 ; [nm] Imaging wavelength
  nd:                3,                    ; padding coefficient for PSF computation
}
{camera_IR_ngs,
  wavelengthInNm:    1064,                 ; [nm] Imaging wavelength
  nd:                3,                    ; padding coefficient for PSF computation
}


; --------------------------- Modal reconstructor ---------------------------
{modalrec,
   intmat_tag:        'auto'                 ; interaction matrix tag
   recmat_tag:        'auto'                 ; reconstruction matrix tag
   nmodes:             19
}
{modalrec_tt, ; to get TT seen by LGS on jitter mirror
  intmat_tag:        'auto'                     ; interaction matrix tag
  recmat_tag:        'auto'                     ; reconstruction matrix tag
  nmodes:             2
}
{modalrec_IR_ngs,
   intmat_tag:        'auto'       ; reconstruction matrix tag
   recmat_tag:        'auto'       ; reconstruction matrix tag
   nmodes:             19
}


; --------------------------- Control ---------------------------
{control,
  delay:             1,                     ; Total temporal delay in time steps
  type:              'INT',                 ; type of control 
  int_gain:         [0.6,0.6,0.4,0.4,replicate(0.3,15)]   ; Integrator gain (for 'INT' control)
}
{control_IR_ngs,
  delay:             1,                     ; Total temporal delay in time steps
  type:              'INT',                 ; type of control
  int_gain:         [0.6,0.6,replicate(0.0,17)]   ; Integrator gain (for 'INT' control)
}


; --------------------------- Atmosphere ---------------------------
{atmo,
   L0:                30,                   ; [m] Outer scale
   heights:           [0,4e3,12e3,20e3] ; Durham 5layer atm
   Cn2:               [0.769,0.104,0.126,0.0] ; Durham 5layer atm
   pixel_phasescreens: 32768L*2                ; size of the phase screen array. Max: 32768L ; Default is 8192L
   fov_in_m      :   8
   infinte_phasescreen : 0B 
}
{seeing,
  constant:         0.9759*0.5/(0.06*4.848)   ; ["] seeing=0.9759*0.5/(r0*4.848)
}
{wind_speed,
  ;constant :        [3.541,84.948,242.644,387.745] ; Durham 5layer atm
  constant :        [3.541,8.399,12.997,5.0] ; Durham 5layer atm
}
{wind_direction,
  constant:         [0.,0.,0.,0.]   ; [degrees] Wind direction value
}



; --------------------------- Pupil stop ---------------------------
{pupil_stop, ; this is used in the calibration
  pupil_mask_tag: 'ALASCA_30cm_120p_down'     
}
{pupil_stop_down,
  pupil_mask_tag: 'ALASCA_30cm_120p_down'     
}
{pupil_stop_up,
  pupil_mask_tag: 'ALASCA_30cm_120p_up'
}


; --------------------------- Modal analysis ---------------------------
{modalanalysis,
   type:              'zernike',             ; modes type
   nmodes:            250,                    ; number of modes
   npixels:           120,                   ; linear dimension of DM phase array
   obsratio:          0.0,                   ; obstruction dimension ratio w.r.t. diameter
}


; --------------------- random jitter of laser spot --------------------

{lgsttres,                ; disturbance to add random Tip/Tilt on LGS SHS
  func_type: 'RANDOM',
  amp: [171,171],     ; total 200 mas RMS on sky --> sqrt(100.^2./2.)*4.848e-9*1.0/4.*1e9
  dm_type: 'zernike',
  nmodes: 2,
  height: 0,
  dm_npixels: 120.,
  dm_obsratio: 0,
  seed: 1
}

{tt_DM, ; zernike DM to remove atmospheric Tip/Tilt from LGS SHS for random jitter of laser spot
  type: 'zernike'
  nmodes: 2
  npixels: 120.,
  diaratio: 1.0, 
  obsratio : 0,
  height: 0
}

{tt_modalanalysis, ; modal analysis to remove atmospheric Tip/Tilt from LGS SHS
  type:             'zernike'
  nmodes:           2,
  npixels:          120.,                        ; linear dimension of DM phase array
  obsratio:         0,                          ; obstruction dimension ratio w.r.t. diameter
}

