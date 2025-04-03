function update_alasca_leo_params, params, objectDir=objectDir, pwfs=pwfs

  rad2arcsec = 3600.*180./!pi
  
  if params.wfs_ngs_source.height lt 38e6 then objectSpeed = 7.6549062e3 $
    else objectSpeed = 3.0e3
  
  ; calculate point ahead angle based on object height and speed
  pointAheadAngleRad = objectSpeed/299792458d
  ;params.wfs_ngs_source.polar_coordinate = [2*pointAheadAngleRad*rad2arcsec,0.]
  ;params.wfs_ngs_source.polar_coordinate = [4.0,0.]
  
  ; consider telescope slewing for LEO case - adapt wind speed and direction based on object speed and direction
;  if params.wfs_ngs_source.height lt 38e6 then begin
;    wind_speed = params.wind_speed.constant
;    wind_direction = params.wind_direction.constant
;    telescope_slewing, (size(params.wind_speed.constant,/dim))[0], params.wfs_ngs_source.height, objectSpeed, $
;    objectDir, pointAheadAngleRad, params.atmo.heights, wind_speed=wind_speed, $
;    wind_direction=wind_direction
;    params.wind_speed.constant = wind_speed
;    params.wind_direction.constant = wind_direction
;  endif
  
  if params.hasKey('extended_object') then extended_object = params.remove('extended_object')

  nmodes = params.modalrec.nmodes
  if ~params.HasKey('calib_disturb_generator') then begin
    if n_elements(lgs_nframes) eq 0 then lgs_nframes = 10L
    if n_elements(ngs_nframes) eq 0 then ngs_nframes = 100L
    if n_elements(seeing4calib) eq 0 then seeing4calib = params.seeing.constant
    aber_name_lgs = '_aber_s'+strtrim(string(seeing4calib,format='(f9.2)'),2)+'asec_fitonly_nm'+strtrim(round(nmodes),2)+'_steps'+strtrim(round(lgs_nframes),2)
    aber_name_ngs = '_aber_s'+strtrim(string(seeing4calib,format='(f9.2)'),2)+'asec_fitonly_nm'+strtrim(round(nmodes),2)+'_steps'+strtrim(round(ngs_nframes),2)
  endif else begin
    aber_name_lgs = ''
    aber_name_ngs = ''
  endelse
  
  
  if pwfs then begin
    ; update PWFS params
    give_me_the_tags_scao, params, params.main.instrument_name, /update_params, new_params=new_params
    params = temporary(new_params)
    if params.modalrec.intmat_tag eq '' then params.modalrec.intmat_tag = params.modalrec.recmat_tag
    if n_elements(extended_object) gt 0 then begin
      if n_elements(extended_object.layerHeight) eq 1 then extTxt = '' else extTxt = '_extObj'
    endif else extTxt = ''
    
    params_ngs = duplicate_params(params)
    params_ngs.wfs_source = params_ngs.wfs_ngs_source
    params_ngs.pyramid = params_ngs.pyramid_IR_ngs
    params_ngs.detector = params_ngs.detector_IR_ngs
    params_ngs.slopec = params_ngs.slopec_IR_ngs
    params_ngs.modalrec = params_ngs.modalrec_IR_ngs

    give_me_the_tags_scao, params_ngs, params.main.instrument_name, /update_params, new_params=new_params_ngs

    params.pyramid_IR_ngs = new_params_ngs.pyramid
    params.wfs_ngs_source = new_params_ngs.wfs_source
    params.slopec_IR_ngs = new_params_ngs.slopec
    params.modalrec_IR_ngs = new_params_ngs.modalrec

    if params.modalrec_IR_ngs.intmat_tag eq '' then params.modalrec_IR_ngs.intmat_tag = params.modalrec_IR_ngs.recmat_tag

    params.modalrec_IR_ngs.intmat_tag = params.modalrec_IR_ngs.intmat_tag+aber_name_ngs
    params.modalrec_IR_ngs.recmat_tag = params.modalrec_IR_ngs.recmat_tag+aber_name_ngs

    if n_elements(extended_object) gt 0 then params.extended_object = extended_object
  endif else begin
    params_sh_lgs = duplicate_params(params)
    params_sh_lgs.sh = params_sh_lgs.sh
    params_sh_lgs.detector = params_sh_lgs.detector
    params_sh_lgs.slopec = params_sh_lgs.slopec

    give_me_the_tags_scao, params_sh_lgs, params.main.instrument_name, /update_params, new_params=new_params_sh_lgs
    params.sh = new_params_sh_lgs.sh
    params.wfs_source = new_params_sh_lgs.wfs_source
    params.slopec = new_params_sh_lgs.slopec
    params.modalrec = new_params_sh_lgs.modalrec

    params_sh_ngs = duplicate_params(params)
    params_sh_ngs.wfs_source = params_sh_ngs.wfs_ngs_source
    params_sh_ngs.sh = params_sh_ngs.sh_ngs
    params_sh_ngs.detector = params_sh_ngs.detector_IR_ngs
    params_sh_ngs.slopec = params_sh_ngs.slopec_IR_ngs
    params_sh_ngs.modalrec = params_sh_ngs.modalrec_IR_ngs

    give_me_the_tags_scao, params_sh_ngs, params.main.instrument_name, /update_params, new_params=new_params_sh
    params.sh_ngs = new_params_sh.sh
    params.wfs_ngs_source = new_params_sh.wfs_source
    params.slopec_IR_ngs = new_params_sh.slopec
    params.modalrec_IR_ngs = new_params_sh.modalrec
  endelse

  return, params

end
