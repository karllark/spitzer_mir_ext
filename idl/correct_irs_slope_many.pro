; program to merge the IRAC data into the DAT_files

pro correct_irs_slope_many,prompt=prompt

readcol,'~/Spitzer/IR_Ext/IRAC/irext_spitzer_ave_phot_subarray_rad3_sky10_20.dat', $
        format='(A15,A6,F8.2,F8.2,F12.4,F12.4,F12.4,F12.4)', $
        name,filter,x,y,flux,flux_unc,sky_flux,sky_flux_unc

uindxs = uniq(name)
n_uindxs = n_elements(uindxs)

dat_path = '~/Python_git/extstar_data/DAT_files/'

tstr = ''
;n_uindxs = 3
for i = 0,(n_uindxs-1) do begin
;    print,name[uindxs[i]]
    newname = strmid(name[uindxs[i]],0,strlen(name[uindxs[i]])-2)
    if (strlen(newname) EQ 7) then newname = repstr(newname,'hd','hd0')
    if (newname EQ 'ngc2024') then newname += '_1'

    corfac_irs = correct_irs_slope(dat_path + newname + '.dat')

    if (not file_test(dat_path + newname + '_irscorfac.dat')) then $
      file_move,dat_path + newname + '.dat',dat_path + newname + '_irscorfac.dat'

    openr,iunit,dat_path + newname + '_irscorfac.dat',/get_lun
    openw,ounit,dat_path + newname + '.dat',/get_lun

    print,newname,'; old corfac_irs'
    tstr = ''
    while (not eof(iunit)) do begin
        readf,iunit,tstr
        if ((strmid(tstr,0,10) EQ 'corfac_irs') AND (strmid(tstr,0,18) NE 'corfac_irs_maxwave')) then begin
            print,tstr
        endif else begin
            printf,ounit,tstr
        endelse
    endwhile

    ; now add the new irs corfac
    print,'new corfac_irs'

    if ((corfac_irs[0] NE 0.0) AND (n_elements(corfac_irs) EQ 3)) then begin
        ostr = 'corfac_irs = ' + strtrim(string(corfac_irs[0],format='(F8.3)'),2)
        print,ostr
        printf,ounit,ostr

        ostr = 'corfac_irs_slope = ' + strtrim(string(corfac_irs[1],format='(F8.3)'),2)
        print,ostr
        printf,ounit,ostr

        ostr = 'corfac_irs_zerowave = ' + strtrim(string(corfac_irs[2],format='(F8.3)'),2)
        print,ostr
        printf,ounit,ostr
    endif else if (corfac_irs[0] NE 0.0) then begin
        ostr = 'corfac_irs = ' + strtrim(string(corfac_irs[0],format='(F8.3)'),2)
        print,ostr
        printf,ounit,ostr
    endif else begin
        print,'no correct_irs measurment possible'
    endelse

    free_lun,iunit
    free_lun,ounit

    if (keyword_set(prompt)) then read,'Continue: ',tstr
endfor

end
