
function correct_irs_slope,dat_file,apply_corfac=apply_corfac,blue_vals=blue_vals

if (keyword_set(apply_corfac)) then no_corfac = 0 else no_corfac = 1

read_star_data,dat_file,star_data,no_corfac=no_corfac

flux4 = get_band_flux(star_data,'IRAC4')
flux15 = get_band_flux(star_data,'IRS15',flux15_unc)
flux24 = get_band_flux(star_data,'MIPS24')

good_to_go = 0

max_good_irs = 0.
if (star_data.irs_exists) then begin
    indxs = where(star_data.irs_spectra.npts GT 0)
    max_good_irs = max(star_data.irs_spectra[indxs].wavelength)
    print,max_good_irs
    if (max_good_irs GT 30.) then good_to_go = 1
endif

if ((good_to_go) AND (flux4 GT 0.) AND (flux24 GT 0.)) then begin
    gindxs = where((star_data.irs_spectra.npts GT 0) AND (star_data.irs_spectra.flux GT 0),n_gindxs)

;    print,star_data.irs_spectra[gindxs].flux

    readcol,'Resp_Curves/080924ch4trans_sub.txt',fwave,fresp
    tflux = interpol(star_data.irs_spectra[gindxs].flux,star_data.irs_spectra[gindxs].wavelength,fwave)
    fflux4 = total(tflux*fresp)/total(fresp)
    kplot,star_data.irs_spectra[gindxs].wavelength,star_data.irs_spectra[gindxs].flux*star_data.irs_spectra[gindxs].wavelength^4,xrange=[5.,41],psym=100, $
          kplot_type='oi'
;    kplot,star_data.irs_spectra[gindxs].wavelength,star_data.irs_spectra[gindxs].flux,xrange=[5.,30],psym=100, $
;          kplot_type='oo'
;        print,flux4,fflux4
    eff_wave4 = total(fwave*tflux*fresp)/total(tflux*fresp)
    print,eff_wave4
    koplot,[eff_wave4],[flux4*eff_wave4^4],psym=1,symsize=2.0
    koplot,[eff_wave4],[fflux4*eff_wave4^4],psym=2,symsize=2.0

    readcol,'Resp_Curves/MIPS24_resp.dat',fwave,fresp
    tflux = interpol(star_data.irs_spectra[gindxs].flux,star_data.irs_spectra[gindxs].wavelength,fwave)
    fflux24 = total(tflux*fresp)/total(fresp)
    eff_wave24 = total(fwave*tflux*fresp)/total(tflux*fresp)
    print,eff_wave24
    koplot,[eff_wave24],[flux24*eff_wave24^4],psym=1,symsize=2.0
    koplot,[eff_wave24],[fflux24*eff_wave24^4],psym=2,symsize=2.0

    print,'8 ratio (phot/spec) = ', flux4/fflux4
    print,'24 ratio (phot/spec) = ', flux24/fflux24

    ; now make the correction
    slope = (flux24/fflux24 - flux4/fflux4)/(eff_wave24 - eff_wave4)
    mod_line = flux4/fflux4 + slope*(star_data.irs_spectra[gindxs].wavelength - eff_wave4)
;    print,mod_line

;    print,fflux4,fflux24

    koplot,star_data.irs_spectra[gindxs].wavelength, $
           mod_line*star_data.irs_spectra[gindxs].flux*(star_data.irs_spectra[gindxs].wavelength^4),psym=100,linestyle=2
;    kplot,star_data.irs_spectra[gindxs].wavelength, $
;           mod_line,psym=100,linestyle=2

    corfac_params = [flux4/fflux4,slope,eff_wave4]

; calculate the IRS blue peakup number
    readcol,'Resp_Curves/bluePUtrans.txt',fwave,fresp
    tflux = interpol(star_data.irs_spectra[gindxs].flux,star_data.irs_spectra[gindxs].wavelength,fwave)
    fflux15 = total(tflux*fresp)/total(fresp)
    eff_wave15 = total(fwave*tflux*fresp)/total(tflux*fresp)
    print,eff_wave15
    koplot,[eff_wave15],[flux15*eff_wave15^4],psym=1,symsize=2.0
    koplot,[eff_wave15],[fflux15*eff_wave15^4],psym=2,symsize=2.0
    if (flux15 GT 0.) then begin
        print,'IRS blue, eff wave = ', eff_wave15
        print,'15 ratio (phot/spec) = ', flux15/fflux15
        blue_vals = [flux15,fflux15,flux15_unc]
    endif else blue_vals = 0.0
    
endif else if ((max_good_irs GT 15.) AND (flux4 GT 0.)) then begin

    gindxs = where((star_data.irs_spectra.npts GT 0) AND (star_data.irs_spectra.flux GT 0),n_gindxs)

;    print,star_data.irs_spectra[gindxs].flux

    readcol,'Resp_Curves/080924ch4trans_sub.txt',fwave,fresp
    tflux = interpol(star_data.irs_spectra[gindxs].flux,star_data.irs_spectra[gindxs].wavelength,fwave)
    fflux4 = total(tflux*fresp)/total(fresp)
    kplot,star_data.irs_spectra[gindxs].wavelength,star_data.irs_spectra[gindxs].flux*star_data.irs_spectra[gindxs].wavelength^4,xrange=[5.,41],psym=100, $
          kplot_type='oi'
;    kplot,star_data.irs_spectra[gindxs].wavelength,star_data.irs_spectra[gindxs].flux,xrange=[5.,30],psym=100, $
;          kplot_type='oo'
;        print,flux4,fflux4
    eff_wave4 = total(fwave*tflux*fresp)/total(tflux*fresp)
    print,eff_wave4
    koplot,[eff_wave4],[flux4*eff_wave4^4],psym=1,symsize=2.0
    koplot,[eff_wave4],[fflux4*eff_wave4^4],psym=2,symsize=2.0

    print,'8 ratio (phot/spec) = ', flux4/fflux4

    ; now make the correction
    corfac = flux4/fflux4

    koplot,star_data.irs_spectra[gindxs].wavelength, $
           corfac*star_data.irs_spectra[gindxs].flux*(star_data.irs_spectra[gindxs].wavelength^4),psym=100,linestyle=2
;    kplot,star_data.irs_spectra[gindxs].wavelength, $
;           mod_line,psym=100,linestyle=2

    corfac_params = [corfac]

endif else begin
    print,'not all the neeed data available for ' + dat_file
    print,'need IRS spectrum, IRAC4, and (optionally) MIPS24 photometry'
    corfac_params = [0.0]
endelse

return,corfac_params

end
