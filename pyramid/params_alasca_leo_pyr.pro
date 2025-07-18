{main,
 root_dir:          '/home/bstadler/passata/ALASCA_LEO/calibration_pyr_Tenerife_elevation30deg_PAA50murad_mod4/', ; Root directory for calibration manager
 instrument_name:   'ALASCA_LEO'               ; this is 'CaNaPy' to be compatible with the CaNaPy calibration of 589nm path
 pixel_pupil:       120,                   ; Linear dimension of pupil phase array
 pixel_pitch:       0.008333,              ; [m] Pitch of the pupil phase array
 pad_size:          2048,                   ; number of pixels for zero padding, if smaller than pixle_pupil no padding applied
 total_time:        1.000d,                ; [s] Total simulation running time
 time_step:         0.0005d,                ; [s] Simulation time step
 zenithAngleInDeg:  60.0,
 precision:         0,
 verbose:           0B
}

; --------------------------- Deformable mirrors ---------------------------
{DM,
    ifunc_tag  : 'CaNaPy_30cm_dm_zonal_ifunc_120pix'        ; modes-to-commands matrix tag
    m2c_tag    : 'CaNaPy_30cm_dm_m2c_68masters_66modes'     ; DM command-to-phase influence matrix tag
    nmodes     : 8
    height     : 0                           ; DM height [m]
}
;{DM_UP,
; type:              'zernike',             ; modes type
; nmodes:            2,                     ; number of modes
; npixels:           120,                   ; linear dimension of DM phase array
; obsratio:          0.0,                   ; obstruction dimension ratio w.r.t. diameter
; height:            0                      ; DM height [m]
;}


; --------------------------- Sources ---------------------------
; LGS for HO correction
{extended_object,
  polar_coordinate: [0.0, 0]
  height: 90000                  ; source altitude in m for atmosphere propagation (!VALUES.F_INFINITY for cylindrical propagation, 90000 for sodium layer cone effect)
  magnitude : 5.44,
  wavelengthInNm: 589
  type : 'FROM_PSF',             ; see obj_type in compute2d method of extended_source object
  sampling_type : 'POLAR',       ; see sampling_type in compute2d method of extended_source object
  size_obj : 0.1                 ; 2D size, for "TOPHAT" type it is the diameter
  multiples_fwhm : 1.0           ; extended object 2D sampling in multiple of lambda/D (DL FWHM)
  show_source: 0b                ; if set display the extended source points
  layerHeight: 90e3 + 1e3*[-13.231, -10.586, -7.942, -5.297, -2.653, -0.8, 2.636, 5.281, 7.925];, 10.570] ; sodium layers altitude in m, a single elements means a 2D source
  intensityProfile: [ 0.0540, 0.1296, 0.1231, 0.1829, 0.2689, 0.1805, 0.0341, 0.0153, 0.0116];, 0]        ; vector of sodium layers intensity (total = 1), a single element means a 2D source
  ttProfile : fltarr(10,2)
  focusHeight: [90e3]             ; sodium layer focus altitude in m
}
{wfs_source,
  polar_coordinate:  [0.0, 0.0],           ; [arcsec, degrees] source polar coordinates
  magnitude:         6.3,                 ; source magnitude. Uses n_phot.pro
  zeroPoint:         390.d-10              ; Uses n_phot.pro
  wavelengthInNm:    589,                  ; [nm] wavelength
  height:            90000.                ; [m] LGS source altitude
}
; Satellite uplink
{wfs_sat_source,
  polar_coordinate:  [0.0, 0.0],           ; [arcsec, degrees] source polar coordinates
  magnitude:         0,                    ; source magnitude. Uses n_phot.pro
  wavelengthInNm:    1055,                 ; [nm] wavelength
  height:            400e3               ; [m] source altitude
}
; Satellite downlink
{wfs_ngs_source,
  polar_coordinate:  [0.0, 0.0],           ; [arcsec, degrees] source polar coordinates
  magnitude:         0,                    ; source magnitude. Uses n_phot.pro
  wavelengthInNm:    1064,                 ; [nm] wavelength
  height:            400e3               ; [m] source altitude
}


; --------------------------- Wavefront sensors ---------------------------
{pyramid,
   pup_diam:          40.                    ; Pupil diameter in subaps.
   pup_dist:          72.                    ; Separation between pupil centers in subaps.
   fov:               8.0                   ; Requested field-of-view [arcsec]
   fov_errinf :       0.1                    ; Maximum error in reducing fov [$$$$ to be checked]
   fov_errsup:        30.0                   ; To avoid FOV error
   mod_amp:           4.0                    ; Modulation radius (in lambda/D units)
   output_resolution: 240                    ; Output sampling [usually corresponding to CCD pixels]
   fft_res:           3.0                    ; pyramid focal-plane PSF sampling in lambda/D units		
   wavelengthInNm:    589                    ; [nm] Pyramid wavelength
}
{slopec,
  pupdata_tag :      'auto',                ; tag of the pyramid WFS pupils
  sn_tag:            'auto'                 ; tag of the slope reference vector
  thr1 = 0.3                                ; threshold no. 1 for pupil computation
  thr2 = 0.3                                ; threshold no. 2 for pupil computation
}
{pyramid_IR_ngs,
   pup_diam:          40.                    ; Pupil diameter in subaps.
   pup_dist:          72.                    ; Separation between pupil centers in subaps.
   fov:               8.0                    ; Requested field-of-view [arcsec]
   fov_errinf :       0.1                    ; Maximum error in reducing fov [$$$$ to be checked]
   fov_errsup:        30.0                   ; To avoid FOV error
   mod_amp:           0.0                    ; Modulation radius (in lambda/D units)
   output_resolution: 240                    ; Output sampling [usually corresponding to CCD pixels]
   fft_res:           3.0                    ; pyramid focal-plane PSF sampling in lambda/D units   
   wavelengthInNm:    1064                    ; [nm] Pyramid wavelength
}
{slopec_IR_ngs,
  pupdata_tag: 'auto',       ; tag of the sub-apertures used by the WFS
  sn_tag     : 'auto',        ; tag of the reference slope vector
  thr1 = 0.3                               ; threshold no. 1 for pupil computation
  thr2 = 0.3                               ; threshold no. 2 for pupil computation
  thr_value    : 0.0                       ; [ph] threshold value used in the slope computation
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
  dt:                0.0005d,                ; [s] Detector integration time
  bandw:             10,                    ; [nm] Sensor bandwidth
  photon_noise:      1b,                    ; activate photon noise
  readout_noise:     1b,                    ; activate readout noise
  readout_level:     0.8,                   ; readout noise in [e-/pix/frame]
  quantum_eff:       1.0,                    ; quantum efficiency * total transmission
  verbose:           0b
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
   nmodes:             8
}
{modalrec_IR_ngs,
   intmat_tag:        'auto'       ; reconstruction matrix tag
   recmat_tag:        'auto'       ; reconstruction matrix tag
   nmodes:             8
}


; --------------------------- Control ---------------------------
{control,
  delay:             1,                     ; Total temporal delay in time steps
  type:              'INT',                 ; type of control 
  int_gain:         [0.,0.,replicate(0.3,6)]   ; Integrator gain (for 'INT' control)
}
{control_IR_ngs,
  delay:             1,                     ; Total temporal delay in time steps
  type:              'INT',                 ; type of control
  ;int_gain:         [0.8,0.8,replicate(0.,17)]    ; Integrator gain (for 'INT' control)
  int_gain:         [0.5,0.5,replicate(0.,6)]    ; optgains 
}


; --------------------------- Atmosphere ---------------------------
{atmo,
   L0:                20,                   ; [m] Outer scale
   heights:           [955.53979, 7300.5816, 12353.543, 16227.363, 21897.079] ; Durham 5layer atm
   Cn2:               [0.85374713, 0.049742997, 0.073054083, 0.021873636, 0.0015821513] ; Durham 5layer atm
   ;heights:           [0,4e3,12e3,20e3] ; TURBO50
   ;Cn2:               [0.769,0.104,0.126,0.0] ; TURBO50
   ;heights:           [124.34, 7.31e3, 1.26e4, 1.65e4, 2.25e4] ; Durham 5layer atm (daytime)
   ;Cn2:               [0.9793, 0.0087, 0.0108, 0.0011, 7.7447e-5] ; Durham 5layer atm (daytime) - weights (total must be eq 1)
   pixel_phasescreens: 2*32768L                ; size of the phase screen array. Max: 32768L ; Default is 8192L
   fov_in_m      :   8
   infinte_phasescreen : 1B 
}
{seeing,
  ;constant:        0.9759*0.5/(0.06*4.848)  ; ["] seeing=0.9759*0.5/(r0*4.848), 0.9759*0.5/(0.06*4.848)
  constant:        2.5
}

{wind_speed,
  ;constant :        [7.62, 12.52, 10.80, 7.56, 10.45] ; Durham 5layer atm
  constant :        [2.9307, 13.5419, 12.6838, 7.2461, 5.7461] ; Durham 5layer atm
  ;constant :        [3.541,8.399,12.997,5.0]; TURBO50
}
{wind_direction,
  constant:          [90.,90.,90.,90.,90.]   ; [degrees] Wind direction value
  ;constant:          [90.,270.,270.,90.,0.]   ; [degrees] Wind direction value
}


; --------------------------- Pupil stop ---------------------------
{pupil_stop, ; this is used in the calibration
  pupil_mask_tag: 'ALASCA_LEO_30cm_120p_down'     
}
{pupil_stop_down,
  pupil_mask_tag: 'ALASCA_LEO_30cm_120p_down'     
}
{pupil_stop_up,
  pupil_mask_tag: 'ALASCA_LEO_30cm_120p_up'
}


; --------------------------- Modal analysis ---------------------------
{modalanalysis,
   type:              'zernike',             ; modes type
   nmodes:            250,                    ; number of modes
   npixels:           120,                   ; linear dimension of DM phase array
   obsratio:          0.0,                   ; obstruction dimension ratio w.r.t. diameter
}
