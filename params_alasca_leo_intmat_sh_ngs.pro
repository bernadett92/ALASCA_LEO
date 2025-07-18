{calib_source,
  polar_coordinate:    [0.0, 0.0],           ; [arcsec, degrees] source polar coordinates
  magnitude:           3,                    ; source magnitude
  wavelengthInNm:      1064                   ; [nm] wavelength
}

{detector,
  readout_noise     : 0b ; readout noise in e-/pix/frame
  photon_noise      : 0b ; activate photon noise
  background_noise  : 0b ; activate sky background noise
  darkcurrent_noise : 0b ; activate dark current noise
  excess_noise      : 0b ; sctivate excess noise of sqrt(2.) characteristic of EMCCDs.
  readout_level     : 0.0 ; readout noise in e-/pix/frame
  background_level  : 0.0 ; sky background noise in e-/pix/frame
  darkcurrent_level : 0.0 ; dark current value in e-/pix/s
}

{calib_disturb_generator,
  func_type          :'PUSHPULL',
  nmodes             : 69,
  influence_function : 'ALASCA_LEO_20cm_dm_zonal_ifunc_120pix'    ; DM command-to-phase influence matrix tag
  m2c_tag            : 'ALASCA_LEO_20cm_dm_m2c_72masters_69modes' ; modes-to-commands matrix tag
  height             : 0,
  amp                : 30.
}
