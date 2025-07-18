;+
; NAME:
;   read_alpao_if
; PURPOSE:
;   Read the ALPAO Influence Functions from disk
; CATEGORY:
;
; CALLING SEQUENCE:
; function read_alpao_if, basename, IFMASK=IFMASK
; INPUTS:
;   basename           filename prefix. Filenames will be generated
;                      appending a number and a '.int' extension to this basename
; OPTIONAL INPUTS:
;   NACTS              Number of actuators. Defaults to 97
;   MASK_ACTS          If set, mask the pupil regions far from the IF
;                      by subtracting this value from the whole pupil
;                      and setting negative values to zero.
;   DISPLAY            If set, display IFs as they are read
;   VERBOSE            If set, enable verbose output
; RETURN VALUE;
;   fltarr(xsize, ysize, nacts): 3d influence function matrix
; KEYWORD OUTPUTS:
;   IFMASK             fltarr(xsize, ysize) 2d pupil mask
; COMMON BLOCKS:
;   None.
; SIDE EFFECTS:
;   None.
; RESTRICTIONS:
;   None
; MODIFICATION HISTORY:
;   Created 10-DEC-2021 by Alfio Puglisi alfio.puglisi@inaf.it
;-

function read_alpao_if, basename, IFMASK=IFMASK, NACTS=NACTS, MASK_ACTS=MASK_ACTS, $
                        VERBOSE=VERBOSE, DISPLAY=DISPLAY

    if not keyword_set(NACTS) then NACTS=97
    if not keyword_set(basename) then $
       basename='ALPAO_IF/CaNaPy_DM_Calib_INT_act.0_'

    for act=0,NACTS-1 do begin

        fname = basename+strtrim(act,2)+'.int'
        if keyword_set(VERBOSE) then print,'Reading '+fname

        img = read_fisba(fname, mask=mask)
        if act eq 0 then begin
            ss = size(img)
            ifmatrix = fltarr(ss[1], ss[2], NACTS)
            ifmask = mask
        endif

        ; actuator masking
        if keyword_set(MASK_ACTS) then begin
           orig_img=img
           img -= MASK_ACTS
           img[where(img lt 0)] = 0
        endif

        ifmatrix[*,*,act] = img
        ifmask *= mask

        if keyword_set(DISPLAY) then tvscl,img
    endfor
    return, ifmatrix
end
