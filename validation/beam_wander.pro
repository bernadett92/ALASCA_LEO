psfs = readfits('/home2/bstadler/ALASCA_LEO/20cmSubPupil_15degElev_25muradPAA/20cm_TURBO50_1s_2kHz_sh_15degElev_25muradPAA_woAO_detailled/psf_sat.fits')
height = 400e3
pitch = 0.008333
wavelength = 1055e-9

s = size(psfs,/dim)
time = s[0]
angle = fltarr(time)
center = [s[1]/2-1,s[2]/2-1]
dxSatellite = (wavelength*height)/(120*pitch);

for t=0,time-1 do begin
  tmp = reform(psfs[t,*,*])
  maxPsf = max(tmp,location)
  location = array_indices(tmp,location) 
  shift_xy = abs(location-center)*dxSatellite ; m
  angle[t] = atan(sqrt(shift_xy[0]^2+shift_xy[1]^2)/height) ; rad
endfor

print, rms(angle)
print, stddev(angle)
end