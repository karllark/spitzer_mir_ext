; program to merge the IRAC data into the DAT_files

pro merge_ir_irac

readcol,'~/Spitzer/IR_Ext/IRAC/irext_spitzer_ave_phot_subarray_rad3_sky10_20.dat', $
        format='(A15,A6,F8.2,F8.2,F12.4,F12.4,F12.4,F12.4)', $
        name,filter,x,y,flux,flux_unc,sky_flux,sky_flux_unc

flux_unc = sqrt(flux_unc^2 + (flux*0.02)^2)

uindxs = uniq(name)
n_uindxs = n_elements(uindxs)

dat_path = '~/Python_git/extstar_data/DAT_files/'

;n_uindxs = 1
for i = 0,(n_uindxs-1) do begin
;    print,name[uindxs[i]]
    newname = strmid(name[uindxs[i]],0,strlen(name[uindxs[i]])-2)
    if (strlen(newname) EQ 7) then newname = repstr(newname,'hd','hd0')
    if (newname EQ 'ngc2024') then newname += '_1'

    if (not file_test(dat_path + newname + '_old.dat')) then $
      file_move,dat_path + newname + '.dat',dat_path + newname + '_old.dat'

    openr,iunit,dat_path + newname + '_old.dat',/get_lun
    openw,ounit,dat_path + newname + '.dat',/get_lun

    print,newname,'; old IRAC'
    tstr = ''
    while (not eof(iunit)) do begin
        readf,iunit,tstr
        if (strmid(tstr,0,4) EQ 'IRAC') then begin
            print,tstr
        endif else begin
            printf,ounit,tstr
        endelse
    endwhile

    ; now add the new IRAC data
    print,'new IRAC'
    gindxs = where(name EQ name[uindxs[i]],n_gindxs)
    for k = 0,(n_gindxs-1) do begin
        ; includes correction for Bohlin, Gordon, et al. (2011, AJ) work
        case filter[gindxs[k]] of
            'IRAC1' : apfac = 1.112*(1.0-0.023)
            'IRAC2' : apfac = 1.113*(1.0-0.019)
            'IRAC3' : apfac = 1.125*(1.0-0.020)
            'IRAC4' : apfac = 1.218*(1.0-0.005)
            else : stop
        endcase
        ostr = filter[gindxs[k]] + ' = ' + string(flux[gindxs[k]]*apfac,format='(E8.2)') + $
               ' +/- ' + string(flux_unc[gindxs[k]]*apfac,format='(E8.2)') + ' mJy'
        print,ostr
        printf,ounit,ostr
    endfor

    free_lun,iunit
    free_lun,ounit
endfor

end
