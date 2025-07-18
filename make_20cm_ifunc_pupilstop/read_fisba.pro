;+
; NAME:
;   read_fisba
; PURPOSE:
;   Read the Fisba interferogram from disk
; CATEGORY:
;
; CALLING SEQUENCE:
; function read_fisba, filename, MASK=MASK
; INPUTS:
;   filename           interferogram filename.
; OPTIONAL INPUTS:
;   None
; RETURN VALUE;
;   fltarr(xsize, ysize): 2d interferogram
; KEYWORD OUTPUTS:
;   MASK             fltarr(xsize, ysize) 2d pupil mask
; COMMON BLOCKS:
;   None.
; SIDE EFFECTS:
;   None.
; RESTRICTIONS:
;   Among the various Fisba output formats, only 16-bit .int is supported
; MODIFICATION HISTORY:
;   Created 10-DEC-2021 by Alfio Puglisi alfio.puglisi@inaf.it
;-
function read_fisba, filename, MASK=MASK

openr, unit, filename, /GET_LUN

; Skip initial comments
line=''
comments=1b
readf, unit, line
while comments do begin
    if strmid(line,0,1) ne '!' then begin
       comments=0b
    endif else begin
       readf, unit, line
    endelse
endwhile

; Read parameter line and extract x/y size
readf, unit, line

fields = strsplit(line, /extract)
xsize = fields[1]
ysize = fields[2]

img = fltarr(1)
while not EOF(unit) do begin
  readf, unit, line
  values = strsplit(line, /extract)
  img = [img, values]
endwhile

close,unit
free_lun, unit

img = img[1:*]
img = reform(img, xsize, ysize)

mask = img*0+1
mask[where(img eq -32768)] = 0
img[where(mask eq 0)]=0
return, img

end
