; program to merge the IRSB data into the DAT_files

pro merge_ir_irsb

readcol,'~/Dust/Ext/IRSB/irext_spitzer_phot_irs_rad12_sky18_25.dat', $
        format='(A15,A6,F8.2,F8.2,F10.2,F10.2,F10.2,F10.2)', $
        name,filter,x,y,flux,flux_unc,sky_flux,sky_flux_unc

uindxs = uniq(name)
n_uindxs = n_elements(uindxs)

;n_uindxs = 1
for i = 0,(n_uindxs-1) do begin
    name[uindxs[i]] = strlowcase(name[uindxs[i]])
    print,name[uindxs[i]]
    newname = strmid(name[uindxs[i]],0,strlen(name[uindxs[i]]))
    if (strlen(newname) EQ 7) then newname = repstr(newname,'hd','hd0')
    if (newname EQ 'ngc2024') then newname += '_1'
    if (newname EQ 'bd+631964') then newname = 'bd+63d1964'

    if (not file_test('DAT_files/' + newname + '_old4.dat')) then $
      file_move,'DAT_files/' + newname + '.dat','DAT_files/' + newname + '_old4.dat'

    openr,iunit,'DAT_files/' + newname + '_old4.dat',/get_lun
    openw,ounit,'DAT_files/' + newname + '.dat',/get_lun

    print,newname,'; old IRS15'
    tstr = ''
    while (not eof(iunit)) do begin
        readf,iunit,tstr
        if (strmid(tstr,0,5) EQ 'IRS15') then begin
            print,tstr
        endif else begin
            printf,ounit,tstr
        endelse
    endwhile

    ; now add the new IRSB data
    print,'new IRS15'
    gindxs = where(name EQ name[uindxs[i]],n_gindxs)
    for k = 0,(n_gindxs-1) do begin
        apfac = 1.077
        col_cor = 1.038
        apfac /= col_cor
        ostr = 'IRS15 = ' + string(flux[gindxs[k]]*apfac,format='(E8.2)') + $
               ' +/- ' + string(flux_unc[gindxs[k]]*apfac,format='(E8.2)') + ' mJy'
        print,ostr
        printf,ounit,ostr
    endfor
    
    free_lun,iunit
    free_lun,ounit
endfor

end
