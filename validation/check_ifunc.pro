set_white_plot

dir = '/home/bstadler/passata/ALASCA_LEO/calibration_sh_Tenerife_elevation15deg_PAA25murad/'
filename = 'ifunc/CaNaPy_30cm_dm_zonal_ifunc_120pix'
ifuncs = readfits(dir+filename+'.fits', exten=1, /silent)
mask = readfits(dir+filename+'.fits', exten=2, /silent)
filename = 'pupilstop/ALASCA_30cm_120p_down'
pupilstop = readfits(dir+filename+'.fits', exten=1, /silent)
print, size(ifuncs,/dim)

dir = '/home/bstadler/passata/ALASCA_LEO/ifunc/ifunc/'
filename = 'ALASCA_LEO_120pix_54nact_zonal_ifs'
ifuncs1 = readfits(dir+filename+'.fits', exten=1, /silent)
mask1 = readfits(dir+filename+'.fits', exten=2, /silent)
dir = '/home/bstadler/passata/ALASCA_LEO/calibration_sh_Tenerife_elevation15deg_PAA25murad_/'
filename = 'pupilstop/ALASCA_20cm_120p_down'
pupilstop1 = readfits(dir+filename+'.fits', exten=1, /silent)
print, size(ifuncs1,/dim)

tmp = mask
tmp(where(mask1 eq 1)) = 0.5
window, 7
image_show, /as, /sh, tmp

tmp = pupilstop
tmp(where(pupilstop1 eq 1)) = 0.5
window, 8
image_show, /as, /sh, tmp

end