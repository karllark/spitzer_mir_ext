; program to merge the IRAC data into the DAT_files

pro compute_ir_corfac_many

readcol,'~/Spitzer/IR_Ext/IRAC/irext_spitzer_ave_phot_subarray_rad3_sky10_20.dat', $
        format='(A15,A6,F8.2,F8.2,F12.4,F12.4,F12.4,F12.4)', $
        name,filter,x,y,flux,flux_unc,sky_flux,sky_flux_unc

uindxs = uniq(name)
n_uindxs = n_elements(uindxs)

;n_uindxs = 3
for i = 0,(n_uindxs-1) do begin
;    print,name[uindxs[i]]
    newname = strmid(name[uindxs[i]],0,strlen(name[uindxs[i]])-2)
    if (strlen(newname) EQ 7) then newname = repstr(newname,'hd','hd0')
    if (newname EQ 'ngc2024') then newname += '_1'

    corfac_irs = compute_ir_corfac('DAT_files/' + newname + '.dat')

    if (not file_test('DAT_files/' + newname + '_old2.dat')) then $
      file_move,'DAT_files/' + newname + '.dat','DAT_files/' + newname + '_old2.dat'

    openr,iunit,'DAT_files/' + newname + '_old2.dat',/get_lun
    openw,ounit,'DAT_files/' + newname + '.dat',/get_lun

    print,newname,'; old corfac_irs'
    tstr = ''
    while (not eof(iunit)) do begin
        readf,iunit,tstr
        if (strmid(tstr,0,4) EQ 'corfac_irs') then begin
            print,tstr
        endif else begin
            printf,ounit,tstr
        endelse
    endwhile

    ; now add the new irs corfac
    print,'new corfac_irs'

    if (corfac_irs NE 1.0) then begin
        ostr = 'corfac_irs = ' + strtrim(string(corfac_irs,format='(F8.3)'),2)
        print,ostr
        printf,ounit,ostr
    endif else begin
        print,'no corfac_irs measurment possible'
    endelse
    
    free_lun,iunit
    free_lun,ounit
endfor

end
