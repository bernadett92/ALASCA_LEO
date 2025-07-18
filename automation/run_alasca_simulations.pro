;+
;NAME :
; run_sdlr_simulation

; Runs a complete (LGS AO + uplink pre-compensation) SDLR simulation with given parameters.
; Performs a postprocessing step for Lumi data exchange.
;
; Input parameter:
; prefix - prefix of the result file and directory
; params_file_night - PASSATA param file
; savedir - result directory
; seed - random seed
;
; Output:
; PASSATA simulation files are stored in savedir/prefix
; Lumi fits file is stored in /home/bstadler/passata/SDLR/results/lumi/prefix.fits
;
; Written by:     B. Stadler; Jan, 2025
;-


; copy from PASSATA lib (otherwise running from command line was not possible)
function read_params_file, filename, verbose=verbose, expand=expand, savedparams=savedparams
  if keyword_set(savedparams) then begin
    if keyword_set(expand) then return, list(savedparams) $
    else return, savedparams
  endif

  dict = dictionary()

  line    = ''
  ss      = ''
  openr,unit,filename,/get_lun
  while not eof(unit) do begin
    readf,unit,line

    ; Change '=' to ':' at keyword definition
    pos=strsplit(line, '[:=]', /regex)
    if n_elements(pos) gt 1 then begin
      strput, line, ':', pos[1]-1
    endif

    ; Remove comments
    if strpos(line, ';') ge 0 then line = strmid(line,0,strpos(line,';'))

    ; Trim lines and skip empty lines
    line = strtrim(line,2)
    if strlen(line) eq 0 then continue

    ; Accept commas, $ continuations and also nothing at all
    last = strmid(line, 0, 1, /reverse)
    if last eq '$' then line=strmid(line, 0, strlen(line)-1)
    if last ne ',' and last ne '$' and last ne '}' then line += ','

    if last ne '}' then ss += line $
    else begin
      namelen     = strpos(ss,',')
      sname       = strmid(ss,1,namelen-1)         ; extract structure name
      scontents   = strmid(ss,namelen+1,strlen(ss)-namelen-2)+'}'           ; extract structure contents without name
      cmd         = 'dict.'+sname+' = dictionary({'+scontents+')'
      if keyword_set(verbose) then print, cmd
      res         = execute(cmd)
      ss          = ''
      if not res then begin
        print, cmd
        message,'Error reading parameters file: '+ filename
      endif
    endelse
  endwhile
  close,unit
  free_lun,unit
  if keyword_set(expand) then return, iterate_dictionary(dict)
  return, dict
end


pro run_alasca_simulations, atmo, params_file, savedir, seed, simul, simul_time, aberration

  ; set paths for calling the procedure from command line
  PASSATA_FOLDER='/home/bstadler/Simulator/PASSATA'
  coyote=filepath(ROOT=[PASSATA_FOLDER], SUB=['IDL','coyote']  , '')
  !PATH=!PATH+':'+EXPAND_PATH('+'+coyote)
  IDLAstro=filepath(ROOT=[PASSATA_FOLDER], SUB=['IDL','IDLAstro']  , '')
  !PATH = !PATH+':'+EXPAND_PATH('+'+IDLAstro)
  IdlTools=filepath(ROOT=[PASSATA_FOLDER], SUB=['IDL','IdlTools'] , '')
  !PATH=!PATH+':'+EXPAND_PATH('+'+IdlTools)
  Optimization=filepath(ROOT=[PASSATA_FOLDER], SUB=['IDL','Optimization'] , '')
  !PATH=!PATH+':'+EXPAND_PATH('+'+Optimization)
  PASSATA=filepath(ROOT=[PASSATA_FOLDER], SUB=['IDL','PASSATA'] , '')
  !PATH=!PATH+':'+EXPAND_PATH('+'+PASSATA)
  PASSATA=filepath(ROOT=[PASSATA_FOLDER], SUB=['mainPASSATA','ALASCA_LEO'] , '')
  !PATH=!PATH+':'+EXPAND_PATH('+'+PASSATA)
  PASSATA=filepath(ROOT=[PASSATA_FOLDER], SUB=['mainPASSATA','ASI_SDLR'] , '')
  !PATH=!PATH+':'+EXPAND_PATH('+'+PASSATA)

  dir = ROUTINE_DIR()

  ; parameter update
  params = read_params_file(dir+params_file)
  params.main.total_time = simul_time
  params.atmo.seed = seed
  params.aberration.constant = [0,0,aberration]
  
  if atmo eq 'TURBO50' then begin
    params.main.root_dir='/home/bstadler/passata/ALASCA_LEO/calibration_20cm_sh_TURBO50_elevation15deg_PAA25murad_1550nm/'
    params.atmo.heights = [0,4e3,12e3,20e3]
    params.atmo.Cn2 = [0.769,0.104,0.126,0.0]
    params.seeing.constant = 0.9759*0.5/(0.06*4.848)
    params.wind_speed.constant = [3.541,8.399,12.997,5.0]
    params.wind_direction.constant = [0.,0.,0.,0.]
  endif
  if atmo eq 'Tenerife' then begin
    params.main.root_dir='/home/bstadler/passata/ALASCA_LEO/calibration_20cm_sh_Tenerife_elevation15deg_PAA25murad_1550nm/'
    params.atmo.heights = [124.34, 7.31e3, 1.26e4, 1.65e4, 2.25e4]
    params.atmo.Cn2 = [0.9793, 0.0087, 0.0108, 0.0011, 7.7447e-5]
    params.seeing.constant = 2.5
    params.wind_speed.constant = [2.9307, 13.5419, 12.6838, 7.2461, 5.7461]
    params.wind_direction.constant = [0.,0.,0.,0.,0.]
  endif
  if atmo ne 'Tenerife' and atmo ne 'TURBO50' then begin
    print, 'Atmosphere not supported!'
  endif
  
  if simul eq 'noAO' then begin
    params.control.int_gain = [0,0,replicate(0,3)]
    params.control_IR_ngs.int_gain = [0,0,replicate(0,3)]
  endif
  if simul eq 'TT' then begin
    params.control.int_gain = [0,0,replicate(0,3)]
    params.control_IR_ngs.int_gain = [0.5,0.5,replicate(0,3)]
  endif
  if simul eq 'SATAO' then begin
    params.control.int_gain = [0,0,replicate(0,3)]
    params.control_IR_ngs.int_gain = [0.5,0.5,replicate(0.3,3)]
  endif
  if simul eq 'LGSAO' then begin
    params.control.int_gain = [0.5,0.5,replicate(0.3,3)]
    params.control_IR_ngs.int_gain = [0.5,0.5,replicate(0,3)]
  endif
  
  if simul ne 'noAO' and simul ne 'TT' and simul ne 'SATAO' and simul ne 'LGSAO' then begin
    print, 'Simulation not supported!'
  endif else begin
     ; output directory name
    tn = atmo + '_' + simul
    tn += '_seed'+strtrim(string(seed,format='(f9.1)'),2)
    tn += '_beamDivergence'+strtrim(string(aberration,format='(f9.1)'),2)
    print, tn

    if ~file_test(savedir+tn+'.fits') then begin
      ; run on GPU
      detect_dlms, /parallel

      alasca_leo_cloop, params, savedir, tn
    endif else print, "File already exists!"
  endelse


end
