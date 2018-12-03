@ra_dec_cnvt

; program to extract the photometry for the Spitzer observations
; of the Webb calibrators

pro get_webbcal_phot_irs,save_png=save_png,small_ap=small_ap

silent = 1
find_peak = 1

readcol,'irs_nl_data.dat',format="(A15,A12)",names,aorids
aorids = 'r'+aorids
n_names = n_elements(names)

if (not keyword_set(small_ap)) then begin
    ext_file_str = '_irs_rad12_sky18_25'
    irs_ap_rad = 12.
    irs_ap_sky = [18.,25.]
endif else begin
    ext_file_str = '_irs_rad6_sky18_25'
    irs_ap_rad = 6.
    irs_ap_sky = [10.,20.]
endelse

openw,unit1,'irext_spitzer_phot'+ext_file_str+'.dat',/get_lun
printf,unit1,'#Spitzer photometry of IR Extinction calibration stars'
printf,unit1,'#Karl Gordon, get_webbcal_phot_subarray.pro, 7 Sep 2011'
printf,unit1,'#'
printf,unit1,'#name, filter, xpos, ypos, flux[mJy], unc[mJy], total sky[mJy], total sky unc[mJy]'

openw,unit2,'irext_spitzer_indiv_phot'+ext_file_str+'.dat',/get_lun
printf,unit2,'#Spitzer photometry of IR Extinction calibration stars'
printf,unit2,'#Derived from individual BCD images'
printf,unit2,'#Karl Gordon, get_webbcal_phot_irs.pro, 7 Sep 2011'
printf,unit2,'#'
printf,unit2,'#name, filenum, filter, xpos, ypos, flux[mJy], unc[mJy], total sky[mJy], total sky unc[mJy] '

; get Spitzer data
;n_names = 1
for k = 0,(n_names-1) do begin
    print,'working on ' + names[k]

    ; see if there is IRS data
    files = file_search('data/*'+names[k]+'*.fits',count=n_files)
    if (n_files GT 0) then begin
        print,'found IRS data n = ' + strtrim(n_files,2)
        fits_read,files[0],image,header
        extast,header,ast_info
        getrot,ast_info,rot,cdelt
        image_scale = 0.5*total(abs(cdelt)*3600.)

        ; get the ra/dec ref from the 1st BCD file
        readcol,'Lis/'+names[k]+'_irspeakb.lis',bfiles,format='(A)'
        fits_read,bfiles[0],pos_image,pos_header
        ra_ref = fxpar(pos_header,'RA_REF')
        dec_ref = fxpar(pos_header,'DEC_REF')

 ;       print,ra_ref,dec_ref
 ;       stop

        ad2xy,ra_ref,dec_ref,ast_info,x_ref,y_ref
        if (keyword_set(save_png)) then make_png = repstr(files[0],'.fits','_phot') else make_png = 0
        if (keyword_set(small_ap)) then begin ; find source
            kgphot,image,x_ref,y_ref,obj_flux,obj_sn,sky_flux,obj_flux_unc, $
                   ap_rad=12.*image_scale,ap_sky=[18.,25.]*image_scale,make_png=make_png, $
                   silent=silent,find_peak=find_peak,image_scale=image_scale
        endif
        kgphot,image,x_ref,y_ref,obj_flux,obj_sn,sky_flux,obj_flux_unc,sky_flux_unc, $
               ap_rad=irs_ap_rad*image_scale,ap_sky=irs_ap_sky*image_scale,make_png=make_png, $
               silent=silent,find_peak=find_peak,image_scale=image_scale
        filter = 'IRSB'
        printf,unit1,names[k],filter,x_ref,y_ref,obj_flux,obj_flux/obj_sn,sky_flux,sky_flux_unc, $
               format='(A15,A6,F8.2,F8.2,E10.2,E10.2,E10.2,E10.2)'

        ; now get the photometry from the individual images
        ifiles = file_search(aorids[k]+'/ch0/bcd/*_bcdb.fits', $
                             count=n_ifiles)

        if (n_ifiles GT 0) then begin
                ; get new ra,dec
            xy2ad,x_ref,y_ref,ast_info,new_ra,new_dec
            for m = 0,(n_ifiles-1) do begin
                fits_read,ifiles[m],iimage,iheader
                iimage_size = size(iimage)
                extast,iheader,iast_info
                getrot,iast_info,irot,icdelt
                iimage_scale = 0.5*total(abs(icdelt)*3600.)
                ad2xy,new_ra,new_dec,iast_info,ix_ref,iy_ref
                if ((ix_ref GE 0) AND (ix_ref LT iimage_size[1]) AND $
                    (iy_ref GE 0) AND (iy_ref LT iimage_size[2])) then begin
                    if (keyword_set(save_png2)) then make_png = repstr(files[0],'.fits','_phot'+strtrim(string(m+1),2)) $
                    else make_png = 0
                    kgphot,iimage,ix_ref,iy_ref,obj_flux,obj_sn,sky_flux,obj_flux_un,sky_flux_unc, $
                           ap_rad=irs_ap_rad*iimage_scale,ap_sky=irs_ap_sky*iimage_scale,make_png=make_png, $
                           silent=silent,find_peak=0,image_scale=iimage_scale
                    filter = 'IRSB'
                    printf,unit2,names[k],m+1,filter,ix_ref,iy_ref,obj_flux,obj_flux/obj_sn,sky_flux,sky_flux_unc, $
                           format='(A15,I4,A6,F8.2,F8.2,E10.2,E10.2,E10.2,E10.2)'
                endif
            endfor
        endif
    endif
endfor
free_lun,unit1
free_lun,unit2

end
