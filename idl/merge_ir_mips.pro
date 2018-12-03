; program to merge the MIPS data into the DAT_files

pro merge_ir_mips

readcol,'~/Dust/Ext/MIPS/irext_mips24_all.dat', $
        format='(A5,I12,A15,F15.3,F12.3,F10.3,F12.3,F12.3,F10.3,F12.3)', $
        tag,aorkey,name,obstime,aflux,asn,afluxmjy,pflux,psn,flux

flux_unc = flux/psn

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

    if (not file_test('DAT_files/' + newname + '_old.dat')) then $
      file_move,'DAT_files/' + newname + '.dat','DAT_files/' + newname + '_old.dat'

    openr,iunit,'DAT_files/' + newname + '_old.dat',/get_lun
    openw,ounit,'DAT_files/' + newname + '.dat',/get_lun

    print,newname,'; old MIPS'
    tstr = ''
    while (not eof(iunit)) do begin
        readf,iunit,tstr
        if (strmid(tstr,0,4) EQ 'MIPS') then begin
            print,tstr
        endif else begin
            printf,ounit,tstr
        endelse
    endwhile

    ; now add the new MIPS data
    print,'new MIPS'
    gindxs = where(name EQ name[uindxs[i]],n_gindxs)
    for k = 0,(n_gindxs-1) do begin
        apfac = 1.0
        ostr = 'MIPS24 = ' + string(flux[gindxs[k]]*apfac,format='(E8.2)') + $
               ' +/- ' + string(flux_unc[gindxs[k]]*apfac,format='(E8.2)') + ' mJy'
        print,ostr
        printf,ounit,ostr
    endfor
    
    free_lun,iunit
    free_lun,ounit
endfor

end
