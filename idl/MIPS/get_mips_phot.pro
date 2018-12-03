
pro get_mips_phot

aper_only = 1
red_phot_aper_filelist,'irext_mips24_extra',aper_only=aper_only, $
  array='24',/no_centroid,fluxconv=4.54e-2, $
  fit_psf='~/Pro/Convolve/PSFs/Pix1/mips_24_3000K.astrom.conv.24.fits', $
  ap_rad=7.0,ap_sky=[20.,32.],ap_corfac=1.61,/flat2

;red_phot_aper_filelist,'irext_mips24_all',aper_only=aper_only, $
;  array='24',/no_centroid,fluxconv=4.54e-2, $
;  fit_psf='~/Pro/Convolve/PSFs/Pix1/mips_24_3000K.astrom.conv.24.fits', $
;  ap_rad=7.0,ap_sky=[20.,32.],ap_corfac=1.61,/flat2


end
