
pro ave_indiv_phot,show_plot=show_plot,small_ap=small_ap,show_histo=show_histo

if (not keyword_set(small_ap)) then begin
    ext_file_str = '_rad10_sky12_20'
endif else begin
    ext_file_str = '_rad3_sky10_20'
endelse

readcol,'irext_spitzer_indiv_phot_irac_subarray'+ext_file_str+'.dat',format='(A15,I4,A6,F8.2,F8.2,F10.2,F10.2,F10.2,F10.2)', $
        name,fid,filter,x,y,flux,flux_unc,sky_flux,sky_flux_unc

openw,unit2,'irext_spitzer_ave_phot_subarray'+ext_file_str+'.dat',/get_lun
printf,unit2,'#Average spitzer photometry of IR Extinction calibration stars'
printf,unit2,'#average of photometry derived from individual BCD images'
printf,unit2,'#Karl Gordon, ave_indiv_phot.pro, 25 May 2011'
printf,unit2,'#'
printf,unit2,'#name, filter, xpos, ypos, flux[mJy], unc[mJy], sky[mJy], sky unc[mJy], s/n, # images'

indxs = where((name EQ 'snap2_1') OR (name EQ 'snap2_2'),n_indxs)
if (n_indxs GT 0) then begin
    name[indxs] = 'snap2'
endif

sindxs = sort(name)

uindxs = uniq(name[sindxs])
uindxs = sindxs[uindxs]
print,name[uindxs]
n_uindxs = n_elements(uindxs)
for i = 0,(n_uindxs-1) do begin
    indxs = where(name EQ name[uindxs[i]],n_indxs)
    sindxs = sort(filter[indxs])
    indxs = indxs[sindxs]
    ufindxs = uniq(filter[indxs])
    n_ufindxs = n_elements(ufindxs)
    for k = 0,(n_ufindxs-1) do begin
        findxs = where(filter[indxs] EQ filter[indxs[ufindxs[k]]])

        ; sigma clip

        cur_indxs = indxs[findxs]
        iter_sigma_clip,flux,cur_indxs,sigma_vals=sigma_vals,sigma_factor=3.,silent=1
        print,sigma_vals[0],median(flux[cur_indxs])

        if (keyword_set(show_plot)) then begin
            kplot,fid[indxs[findxs]],flux[indxs[findxs]],psym=1, $
                  title=name[uindxs[i]]+ ' / ' + filter[indxs[ufindxs[k]]]
        
            koplot,fid[cur_indxs],flux[cur_indxs],psym=2
            koplot,[-1000,1000],[1.,1.]*sigma_vals[0]

            ans = ''
            read,'Continue: ',ans
        endif

        if (keyword_set(show_histo)) then begin
            min_val = min(flux[cur_indxs],max=max_val)
            n_cur_indxs = n_elements(cur_indxs)
            nbins = max([5,fix(n_cur_indxs/5.)])
            histo_vals = histogram(flux[cur_indxs],min=min_val,max=max_val,nbins=nbins)
            bin_vals = min_val + (max_val-min_val)*(findgen(nbins)+0.5)/float(nbins)

            kplot,bin_vals,histo_vals, $
                  title=name[uindxs[i]]+ ' / ' + filter[indxs[ufindxs[k]]]
            koplot,[1.,1.]*sigma_vals[0],[-1e6,1e6],linestyle=1

            ans = ''
            read,'Continue: ',ans
        endif

;        print,'ave,dev = ',sigma_vals[0],sigma_vals[1]
        x_ref = mean(x[cur_indxs])
        y_ref = mean(y[cur_indxs])
        ave_sky = mean(sky_flux[cur_indxs])
        ave_sky_unc = mean(sky_flux_unc[cur_indxs])

        printf,unit2,name[uindxs[i]],filter[indxs[ufindxs[k]]],x_ref,y_ref,sigma_vals[0], $
              sigma_vals[1],ave_sky,ave_sky_unc,sigma_vals[0]/sigma_vals[1],sigma_vals[3], $
          format='(A15,A6,F8.2,F8.2,E12.4,E12.4,E12.4,E12.4,F8.2,I6)'
    endfor

endfor
free_lun,unit2

end
