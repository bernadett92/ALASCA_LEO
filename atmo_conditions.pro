lamda = 500e-9

; night time
cn2=[0.85374713, 0.049742997, 0.073054083, 0.021873636, 0.0015821513]
h = [955.53979, 7300.5816, 12353.543, 16227.363, 21897.079]

seeing=0.7d
r0= 0.98*lamda/(seeing*!PI/(180.*3600.))
iso_angle_arcsec = iso_planatic_angle(cn2, h, seeing, 500)
iso_angle_rad = iso_angle_arcsec*!PI/(180.*3600)
print, r0
print, iso_angle_rad

seeing=1.0
r0= 0.98*lamda/(seeing*!PI/(180.*3600.))
iso_angle_arcsec = iso_planatic_angle(cn2, h, seeing, 500)
iso_angle_rad = iso_angle_arcsec*!PI/(180.*3600)
print, r0
print, iso_angle_rad

r0=0.1
seeing = (0.98*lamda/r0)*180.*3600./!PI
iso_angle_arcsec = iso_planatic_angle(cn2, h, seeing, 500)
iso_angle_rad = iso_angle_arcsec*!PI/(180.*3600)
print, r0
print, iso_angle_rad

r0=0.05
seeing = (0.98*lamda/r0)*180.*3600./!PI
iso_angle_arcsec = iso_planatic_angle(cn2, h, seeing, 500)
iso_angle_rad = iso_angle_arcsec*!PI/(180.*3600)
print, r0
print, iso_angle_rad

; day time
cn2=[0.9793, 0.0087, 0.0108, 0.0011, 7.7447e-5]
h=[124.34, 7.31e3, 1.26e4, 1.65e4, 2.25e4]

seeing=2.5
r0= 0.98*lamda/(seeing*!PI/(180.*3600.))
iso_angle_arcsec = iso_planatic_angle(cn2, h, seeing, 500)
iso_angle_rad = iso_angle_arcsec*!PI/(180.*3600)
print, r0
print, iso_angle_rad

r0=0.025
seeing = (0.98*lamda/r0)*180.*3600./!PI
iso_angle_arcsec = iso_planatic_angle(cn2, h, seeing, 500)
iso_angle_rad = iso_angle_arcsec*!PI/(180.*3600)
print, r0
print, iso_angle_rad

end