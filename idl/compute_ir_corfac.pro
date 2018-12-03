
function compute_ir_corfac,dat_file

read_star_data,dat_file,star_data,/no_corfac

if (star_data.irs_exists) then begin
    flux4 = get_band_flux(star_data,'IRAC4')
    if (flux4 GT 0.) then begin
        readcol,'Resp_Curves/080924ch4trans_sub.txt',fwave,fresp
;        tfresp = interpol(fresp,fwave,star_data.irs_spectra.wavelength)
;        k_interpolate,fresp,fwave,star_data.irs_spectra.wavelength,tfresp
        tflux = interpol(star_data.irs_spectra.flux,star_data.irs_spectra.wavelength,fwave)
        fflux4 = total(tflux*fresp)/total(fresp)
        kplot,star_data.irs_spectra.wavelength,star_data.irs_spectra.flux,xrange=[5.,10.],psym=100
;        print,flux4,fflux4
        eff_wave = total(fwave*tflux*fresp)/total(tflux*fresp)
        print,eff_wave
        koplot,[eff_wave],[flux4],psym=1,symsize=2.0
        koplot,[eff_wave],[fflux4],psym=2,symsize=2.0
        corfac = flux4/fflux4
    endif else corfac = 1.0
endif else corfac = 1.0

return,corfac

end
