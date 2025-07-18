set_white_plot
display = 1B

savedir = '/home/bstadler/passata/ALASCA_LEO/pupilstop/'

file_full = 'ALASCA_20cm_120p_down.fits'
file_up = 'ALASCA_20cm_120p_up.fits'

npixel = 120
pixel_pitch = 0.008333
diameter = npixel*pixel_pitch
diameter_launch = 0.2
beam_waist = 0.09
obs = 0.29
pixel_shift = floor(npixel*0.5-(npixel-obs*npixel)*0.25+2)

FWHM = 1.18*beam_waist

round_mask_200mm = shift(float(make_mask(npixel,obs=0,dia=diameter_launch/diameter)),pixel_shift,0)
spot_size_norm_pix = FWHM/pixel_pitch
mask_200mm = gaussian_function([spot_size_norm_pix,spot_size_norm_pix]/(2*sqrt(2*alog(2))), WIDTH=npixel, MAXIMUM=1)
mask_200mm = shift(mask_200mm, pixel_shift,0)
mask_200mm *= round_mask_200mm
mask_200mm *= 1/max(mask_200mm)

set_white_plot
;window, 2
;plot, (findgen(npixel)/npixel-0.5), mask_200mm[npixel/2,*], xra=[-0.151,0.151], /xst, thick=2
;oplot, [-1,1], replicate(0.5,2), line=2
;oplot, [-1,1], replicate(1/exp(2),2), line=2
;
;window, 3
;plot, mask_200mm[npixel/2.,*]
;oplot, [0,npixel], [0.5,0.5], linest=2
;oplot, replicate(npixel/2.-0.5-spot_size_norm_pix/2,2),[0,1],linest=2
;oplot, replicate(npixel/2.-0.5+spot_size_norm_pix/2,2),[0,1],linest=2

window, 1
image_show, /as, /sh, round_mask_200mm

window, 2
image_show, /as, /sh, mask_200mm


stop

obj = obj_new('pupilstop', npixel, npixel, pixel_pitch, 0, input_mask=float(round_mask_200mm))

sxaddpar, hdr, 'PUP_STOP', 'ALASCA_LEO'
obj.save, savedir+file_full, hdr

obj = obj_new('pupilstop', npixel, npixel, pixel_pitch, 0, input_mask=float(mask_200mm))

sxaddpar, hdr, 'PUP_STOP', 'ALASCA_LEOUp'
obj.save, savedir+file_up, hdr

end