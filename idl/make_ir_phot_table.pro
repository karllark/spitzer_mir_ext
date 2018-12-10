
pro make_ir_phot_table

readcol,'ir_standards_culled.dat',irs_standards_in,format='(A)'
sindxs = sort(irs_standards_in)
irs_standards_in = irs_standards_in[sindxs]

dat_path = '~/Python_git/extstar_data/DAT_files/'
files = dat_path + irs_standards_in + '.dat'
n_files = n_elements(files)

openw,unit1,'ir_ext_phot.tex',/get_lun
printf,unit1,'\begin{deluxetable*}{lcccccc}'
printf,unit1,'\tablewidth{0pt}'
printf,unit1,'\tablecaption{Spitzer Photometry\label{tab_phot}}'
printf,unit1,'\tablehead{\colhead{Name} & \colhead{IRAC1} & \colhead{IRAC2} & \colhead{IRAC3} & \colhead{IRAC4} & \colhead{IRSB} & \colhead{MIPS24} }'
printf,unit1,'\startdata'
printf,unit1,'\multicolumn{7}{c}{Comparison Stars} \\ \hline'
for i = 0,(n_files-1) do begin
    read_star_data,files[i],star_data
    ostr = strupcase(irs_standards_in[i])

    mag = get_band_mag(star_data,'IRAC1',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'IRAC2',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'IRAC3',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'IRAC4',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'IRS15',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'MIPS24',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    ostr += ' \\'
    printf,unit1,ostr
endfor

readcol,'ir_reddened_av_sort_culled.dat',irs_reddened_in,format='(A)'
sindxs = sort(irs_reddened_in)
irs_reddened_in = irs_reddened_in[sindxs]

files = 'DAT_files/' + irs_reddened_in + '.dat'
n_files = n_elements(files)

printf,unit1,'\multicolumn{7}{c}{Reddened Stars} \\ \hline'
for i = 0,(n_files-1) do begin
    read_star_data,files[i],star_data
    ostr = strupcase(irs_reddened_in[i])

    mag = get_band_mag(star_data,'IRAC1',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'IRAC2',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'IRAC3',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'IRAC4',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'IRS15',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'MIPS24',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    ostr += ' \\'
    printf,unit1,ostr
endfor

printf,unit1,'\enddata'
printf,unit1,'\end{deluxetable*}'
free_lun,unit1

print, 'done'

end
