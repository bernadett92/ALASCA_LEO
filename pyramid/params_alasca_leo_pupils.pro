{calib_source,
 polar_coordinate:    [0.0, 0.0],           ; [arcsec, degrees] source polar coordinates
 magnitude:           2,                    ; source magnitude
 wavelengthInNm:      589                   ; [nm] wavelength
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

{pyramid,
  mod_amp: 4
}
