{main,
 root_dir:          '/home/bstadler/passata/ALASCA_LEO/calibration_20cm_sh_TURBO50_elevation15deg_PAA25murad_1550nm/', ; Root directory for calibration manager
 instrument_name:   'ALASCA_LEO'               ; this is 'CaNaPy' to be compatible with the CaNaPy calibration of 589nm path
 pixel_pupil:       120,                   ; Linear dimension of pupil phase array
 pixel_pitch:       0.008333,              ; [m] Pitch of the pupil phase array
 total_time:        1.00d,                ; [s] Total simulation running time
 time_step:         0.0005d,                ; [s] Simulation time step
 zenithAngleInDeg:  75.0,
 precision:         0,
 verbose:           0B
}

; --------------------------- Deformable mirrors ---------------------------
{DM,
    ifunc_tag  : 'ALASCA_LEO_20cm_dm_zonal_ifunc_120pix'        ; modes-to-commands matrix tag
    m2c_tag    : 'ALASCA_LEO_20cm_dm_m2c_72masters_69modes'     ; DM command-to-phase influence matrix tag
    nmodes     : 69
    height     : 0                           ; DM height [m]
}

; --------------------------- Sources ---------------------------
; Laser ranging uplink
{wfs_sat_source,
  polar_coordinate:  [0.0, 0.0],              ; [arcsec, degrees] source polar coordinates
  magnitude:         0,                       ; source magnitude. Uses n_phot.pro
  wavelengthInNm:    1550,                    ; [nm] wavelength
  height:            400e3                    ; [m] source altitude
}
; Laser ranging downlink
{wfs_ngs_source,
  ;polar_coordinate:  [5.2667691952191076, 180.],  ; [arcsec, degrees] source polar coordinates
  polar_coordinate:  [7.5,0.],  ; [arcsec, degrees] source polar coordinates
  magnitude:         0,                       ; source magnitude. Uses n_phot.pro
  wavelengthInNm:    1064,                    ; [nm] wavelength
  height:            400e3                    ; [m] source altitude
}
; LGS for HO correction
{wfs_source,
  polar_coordinate:  [0., 0.0],              ; [arcsec, degrees] source polar coordinates
  magnitude:         6.3,                     ; source magnitude. Uses n_phot.pro
  zeroPoint:         390.d-10                 ; Uses n_phot.pro
  wavelengthInNm:    589,                     ; [nm] wavelength
  height:            90000.                   ; [m] LGS source altitude
}

; --------------------------- Cameras ---------------------------
{camera_IR_sat,
  wavelengthInNm:    1064,                 ; [nm] Imaging wavelength
  nd:                3,                    ; padding coefficient for PSF computation
}
{camera_IR_ngs,
  wavelengthInNm:    1064,                 ; [nm] Imaging wavelength
  nd:                3,                    ; padding coefficient for PSF computation
}


{control,
  delay:             1,                     ; Total temporal delay in time steps
  type:              'INT',                 ; type of control
  int_gain:         [0.5,0.5,replicate(0.3,67)]   ; Integrator gain (for 'INT' control)
}

; --------------------------- Atmosphere ---------------------------
{atmo,
  L0:                20,                   ; [m] Outer scale
  ;heights:           [124.34, 7.31e3, 1.26e4, 1.65e4, 2.25e4] ; Durham 5layer atm (daytime)
  ;Cn2:               [0.9793, 0.0087, 0.0108, 0.0011, 7.7447e-5] ; Durham 5layer atm (daytime) - weights (total must be eq 1)
  heights:           [0,4e3,12e3,20e3] ; TURBO50
  Cn2:               [0.769,0.104,0.126,0.0] ; TURBO50
   pixel_phasescreens: 2*32768L                ; size of the phase screen array. Max: 32768L ; Default is 8192L
   fov_in_m      :   8
   infinte_phasescreen : 1B 
}
{seeing,
  constant:        0.9759*0.5/(0.06*4.848)  ; ["] seeing=0.9759*0.5/(r0*4.848), 0.9759*0.5/(0.06*4.848)
  ;constant:        2.5
}

{wind_speed,
  ;constant :        [5.3102274,153.43532,253.81335,323.01099,436.33456] ; Durham 5layer atm
  ;constant :        [2.3795276,139.89342,241.12955,315.76489,430.58847]
  ;constant :        [2.9307, 13.5419, 12.6838, 7.2461, 5.7461] ; Durham 5layer atm
  constant :        [3.541,84.948,242.644,387.745]  ; TURBO50
}
{wind_direction,
  constant:          [0.,0.,0.,0.]   ; [degrees] Wind direction value
  ;constant:          [90.,270.,270.,90.,0.]   ; [degrees] Wind direction value
}


; --------------------------- Pupil stop ---------------------------
{pupil_stop, ; this is used in the calibration
  pupil_mask_tag: 'ALASCA_20cm_120p_down'
}
{pupil_stop_down,
  pupil_mask_tag: 'ALASCA_20cm_120p_down'
}
{pupil_stop_up,
  pupil_mask_tag: 'ALASCA_20cm_120p_up'
}


; --------------------------- Modal analysis ---------------------------
{modalanalysis,
  phase2modes_tag: 'ALASCA_LEO_20cm_dm_ifunc_120pix_72masters_69modes_inv'
}