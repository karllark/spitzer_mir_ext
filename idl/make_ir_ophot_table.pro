
pro make_ir_ophot_table

readcol,'ir_standards_culled.dat',irs_standards_in,format='(A)'
sindxs = sort(irs_standards_in)
irs_standards_in = irs_standards_in[sindxs]

files = 'DAT_files/' + irs_standards_in + '.dat'
n_files = n_elements(files)

openw,unit1,'ir_ext_ophot.tex',/get_lun
;printf,unit1,'\begin{deluxetable*}{lcccccccc}'
printf,unit1,'\begin{deluxetable*}{lcccccc}'
printf,unit1,'\tablewidth{0pt}'
printf,unit1,'\tablecaption{Optical/NIR Photometry\label{tab_ophot}}'
;printf,unit1,'\tablehead{\colhead{Name} & \colhead{U} & \colhead{B} & \colhead{V} & \colhead{R} & \colhead{I} & \colhead{J} & \colhead{H} & \colhead{K} }'
printf,unit1,'\tablehead{\colhead{Name} & \colhead{U} & \colhead{B} & \colhead{V} & \colhead{J} & \colhead{H} & \colhead{K} }'
printf,unit1,'\startdata'
;printf,unit1,'\multicolumn{9}{c}{Comparison Stars} \\ \hline'
printf,unit1,'\multicolumn{7}{c}{Comparison Stars} \\ \hline'
for i = 0,(n_files-1) do begin
    read_star_data,files[i],star_data
    ostr = strupcase(irs_standards_in[i])

    mag = get_band_mag(star_data,'U',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'B',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'V',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

;    mag = get_band_mag(star_data,'R',mag_unc)
;    if (mag GT -1.0) then $
;      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
;    else ostr += ' & \nodata'

;    mag = get_band_mag(star_data,'I',mag_unc)
;    if (mag GT -1.0) then $
;      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
;    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'J',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'H',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'K',mag_unc)
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

    mag = get_band_mag(star_data,'U',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'B',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'V',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

;    mag = get_band_mag(star_data,'R',mag_unc)
;    if (mag GT -1.0) then $
;      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
;    else ostr += ' & \nodata'

;    mag = get_band_mag(star_data,'I',mag_unc)
;    if (mag GT -1.0) then $
;      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
;    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'J',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'H',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    mag = get_band_mag(star_data,'K',mag_unc)
    if (mag GT -1.0) then $
      ostr += ' & $' + string(mag,format='(F6.3)') + ' \pm ' + string(mag_unc,format='(F6.3)') + '$' $
    else ostr += ' & \nodata'

    ostr += ' \\'
    printf,unit1,ostr
endfor

printf,unit1,'\enddata'
printf,unit1,'\end{deluxetable*}'
free_lun,unit1

end
