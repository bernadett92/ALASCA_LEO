;+
; NAME:
;   make_alpao_base
; PURPOSE:
;   Generate a modal basis from an ALPAO influence function matrix
; CATEGORY:
;
; CALLING SEQUENCE:
; function make_alpao_base, ifunc_matrix, ifunc_mask, $
;                           AMP=AMP, M2C=M2C, DOUBLE=DOUBLE
; INPUTS:
;   ifunc_matrix:     3d influence function matrix [x,y,acts]
;   ifunc_mask:       2d pupil mask [x,y]
;   amp:              actuation amplitude for each actuator,
;                     or a scalar value identical to all actuators
; OPTIONAL INPUTS:
;   DOUBLE:           perform calculations in double precision
;   VERBOSE:          enable verbose output
;   BAD_ACTS:         vector of actuators to be removed.
;
; RETURN VALUE;
;   fltarr(xsize, ysize,nmodes): 3d modal basis
; KEYWORD OUTPUTS:
;   M2C             fltarr(nmodes, nacts) modes-to-command matrix
; COMMON BLOCKS:
;   None.
; SIDE EFFECTS:
;   None.
; RESTRICTIONS:
;   None
; MODIFICATION HISTORY:
;   Created 10-DEC-2021 by Alfio Puglisi alfio.puglisi@inaf.it
;-

function make_alpao_base, ifunc_matrix, ifunc_mask, amp, $
                          M2C=M2C, DOUBLE=DOUBLE, VERBOSE=VERBOSE, $
                          DISPLAY=DISPLAY, BAD_ACTS=BAD_ACTS

 if not keyword_set(DOUBLE) then DOUBLE=0

 ss = size(ifunc_matrix, /dim)
 nacts = ss[2]
 original_nacts = nacts

 if n_elements(amp) eq 1 then amp =replicate(amp, nacts)

 ; Make sure ifunc images are square
 if ss[0] ne ss[1] then begin
     minsize = min(ss[0:1])
     ifunc = ifunc_matrix[0:minsize-1, 0:minsize-1, *]
     mask = ifunc_mask[0:minsize-1, 0:minsize-1]
     if keyword_set(VERBOSE) then $
         print,'IF truncated from '+strtrim(ss[0],2)+'x'+strtrim(ss[1],2)+' to '+strtrim(minsize,2)+'x'+strtrim(minsize,2)
 endif else begin
     ifunc = ifunc_matrix
     mask = ifunc_mask
 endelse

 ; Amplitude normalization
 for act=0, nacts-1 do ifunc[*,*,act] /= amp[act]

 ; Bad acts removal
 if keyword_set(BAD_ACTS) then begin
    idx = lonarr(nacts)+1
    idx[BAD_ACTS]=0
    use_acts = where(idx)
    ifunc = ifunc[*,*,use_acts]
    nacts = n_elements(use_acts)
 endif


 ; Actuator selection
 ; Do not use, because it messes up the M2C actuator ordering
 if 0 then begin
   mm = fltarr(nacts)
   for i=0,nacts-1 do mm[i] = max(ifunc[*,*,i])-min(ifunc[*,*,i])
   plot,mm
   thr = 30000
   goodacts = where(mm gt thr)
   help,goodacts
   ifunc = ifunc[*,*,goodacts]
   use_acts = use_acts[goodacts]
   nacts = n_elements(use_acts)
 endif
 
 ; Reshape to 2d since make_modal_base_from_ifs_fft.pro
 ; has been tested only with 2d IFs
 ss = size(ifunc, /dim)
 good = where(mask)
 ifunc2d = fltarr(ss[2], n_elements(good))
 for act=0,ss[2]-1 do begin
    img = reform(ifunc[*,*,act], ss[0]*ss[1])
    ifunc2d[act, *] = img[good]
 endfor

 ; Some default parameters...
 diam = 4.0
 r0 = 0.5
 L0 = 40.0
 ZernModes=2
 oversampling=2
 IFmaxConditionNumber = 2e3


 klbasis = make_modal_base_from_ifs_fft( mask, diam, ifunc2d, $
                        r0, L0, ZernModes=ZernModes, oversampling=oversampling, M2C=M2C, $
                        IFmaxConditionNumber=IFmaxConditionNumber, VERBOSE=VERBOSE, DOUBLE=DOUBLE)

 ; Adjust M2C for the removed bad actuators
 nmodes = n_elements(m2c[*,0])
 nacts  = n_elements(m2c[0,*])

 m2c_full = fltarr(nmodes, original_nacts)
 m2c_full[*,use_acts] = m2c
 m2c = m2c_full

 ; Convert KL basis back to 3d
 klsize = size(klbasis,/dim)
 nacts = klsize[0]
 klbasis3d = fltarr(nacts,ss[0],ss[1])
 for i=0,nacts-1 do begin
     img2d = fltarr(ss[0],ss[1])
     img2d[where(mask)] = klbasis[i,*]
     klbasis3d[i,*,*] = img2d
 endfor
 klbasis = klbasis3d

if keyword_set(DISPLAY) then begin
   window,0,xsize=800,ysize=800
   for i=0,nacts-1 do tv, bytscl( congrid(reform(klbasis[i,*,*]),80,80,150,/inter),min=-2, max=2),i
endif

return, klbasis
end
