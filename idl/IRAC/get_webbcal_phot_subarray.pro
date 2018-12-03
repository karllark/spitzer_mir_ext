@ra_dec_cnvt

; program to extract the photometry for the Spitzer observations
; of the Webb calibrators

pro get_webbcal_phot_subarray,save_png=save_png,small_ap=small_ap,test_centering=test_centering,apcor=apcor, $
                              small_sky=small_sky

silent = 1
find_peak = 1

readcol,'irac_nl_data.dat',format="(A12,A12)",names,aorid

print,names
print,aorid

;if (keyword_set(apcor)) then names = ['1802271_' + strtrim(indgen(19)+1,2)]
;if (keyword_set(apcor)) then names = ['1812095_' + strtrim(indgen(117)+1,2)]
;if (keyword_set(apcor)) then names = ['hd165459_' + strtrim(indgen(141)+1,2)]

n_names = n_elements(names)

; bad regions to mask before photometry
bad_reg = intarr(n_names)
bad_reg_file = strarr(n_names)
;bad_reg[0] = 1
;bad_reg_file[0] = 'g191b2b_1_bad_ds9.reg'
;bad_reg[1] = 1
;bad_reg_file[1] = 'g191b2b_2_bad_ds9.reg'
;n_names = 1

; IRAC correction coefficients (PSF)
a = fltarr(6,4)
a[*,0] = [1.0114,-3.536E-6,-6.826E-5,-1.618E-8,1.215E-6,1.049E-6]
a[*,1] = [1.0138,8.401E-5,3.345E-7,1.885E-7,1.438E-6,1.337E-6]
a[*,2] = [1.0055,-3.870E-4,4.600E-5,1.956E-7,2.078E-6,9.970E-7]
a[*,3] = [1.0054,2.332E-4,-8.234E-5,-1.881E-7,6.520E-7,9.415E-7]

; IRAC correction coefficients (pixel phase)
irac_phase_a = [0.0535,0.0309]

if (keyword_set(small_sky)) then begin
    irac_ap_rad = 3.
    irac_ap_sky = [3.,7.]
    ext_file_str = '_irac_subarray_rad3_sky3_7'
endif else if (keyword_set(small_ap)) then begin
    irac_ap_rad = 3.
    irac_ap_sky = [10.,20.]
    ext_file_str = '_irac_subarray_rad3_sky10_20'
endif else begin
    irac_ap_rad = 10.
    irac_ap_sky = [12.,20.]
    ext_file_str = '_irac_subarray_rad10_sky12_20'
endelse
if (keyword_set(apcor)) then ext_file_str = '_apcor' + ext_file_str

; read in IRAC AOR ids
;readcol,'irac_aorids.txt',xnames,xid1,xid2,xid3,xid4,format='(A,A,A,A,A)'
;n_xnames = n_elements(xnames)
;xid = strarr(4,n_xnames)
;xid[0,*] = xid1
;xid[1,*] = xid2
;xid[2,*] = xid3
;xid[3,*] = xid4

;irac_aorids = ['r28510976','r28511232','r28812544','r28812800','r28957696', $
;               'r7658240','r7655680','r7654912','r7656960','r11948800','r25130240']
;irac_indx = [1,5,2,3,0,4,0,11,10,6,7,8,9,0] - 1
;irs_aorids = ['r28647424','r28647680','r28647936','r28648192','r28648448','r28648704','r14912256']
;irs_indx = [0,0,0,0,7,6,0,0,0,1,2,3,4,5,0,0] - 1

;irac_indx = replicate(-1,n_names)
;irs_indx = replicate(-1,n_names)

; get coordinates

;readcol,'coordinates.txt',cor_name,ra_str,dec_str,format='(A12,A12,A12)'
;cor_name = strlowcase(cor_name)

openw,unit1,'irext_spitzer_phot'+ext_file_str+'.dat',/get_lun
printf,unit1,'#Spitzer photometry of IR Extinction calibration stars'
printf,unit1,'#Karl Gordon, get_webbcal_phot_subarray.pro, 25 May 2011'
printf,unit1,'#'
printf,unit1,'#name, filter, xpos, ypos, flux[mJy], unc[mJy], total sky[mJy], total sky unc[mJy]'

openw,unit2,'irext_spitzer_indiv_phot'+ext_file_str+'.dat',/get_lun
printf,unit2,'#Spitzer photometry of IR Extinction calibration stars'
printf,unit2,'#Derived from individual BCD images'
printf,unit2,'#Karl Gordon, get_webbcal_phot_subarray.pro, 25 May 2011'
printf,unit2,'#'
printf,unit2,'#name, filenum, filter, xpos, ypos, flux[mJy], unc[mJy], total sky[mJy], total sky unc[mJy] '

; get Spitzer data
;n_names = 1
for k = 0,(n_names-1) do begin
;for k = (n_names-2),(n_names-1) do begin
;for k = 6,6 do begin
    print,'working on ' + names[k]

    ; get coordinates
;    us_pos = strpos(names[k],'_')
;    if (us_pos GT 0) then begin
;        tcor_name = strmid(names[k],0,us_pos)
;    endif else begin
;        tcor_name = names[k]
;    endelse
;    indxs = where(tcor_name EQ cor_name,n_indxs)
;    if (n_indxs EQ 0) then begin
;        print,'no cor match for ' + tcor_name
;        stop
;    endif
;    ra_ref = 15.*ra_str2hr(ra_str[indxs[0]])
;    dec_ref = dec_str2deg(dec_str[indxs[0]])

    ; see if there is IRAC data
    files = file_search('Lis/*'+names[k]+'*_bcd_exptime0.0?.lis',count=n_files)
;    files = file_search('data/*'+names[k]+'*.fits',count=n_files)
    if (n_files GT 0) then begin
        print,' '
        print,'found IRAC data n = '+strtrim(n_files,2)
        for i = 0,(n_files-1) do begin
            readcol,files[i],bfiles,format='(A)'
            ; read the 1st stack and make a coadded image
            n_bfiles = n_elements(bfiles)

            for z = 0,(n_bfiles-1) do begin
                fits_read,bfiles[z],cube,header
                size_cube = size(cube)
                image = median(cube,dimension=3)

                ra_ref = fxpar(header,'RA_REF')
                dec_ref = fxpar(header,'DEC_REF')

                if (z EQ 0) then begin
                    output_image = image
                    output_header = header
                endif else begin
                    output_image += image
                endelse

                save_image = image
;            if (bad_reg[k]) then begin
;                get_exclude_regions,'data/'+bad_reg_file[k],exclude_regions,n_regions,silent=silent
;                nan_dce_body,image,header,exclude_regions=exclude_regions,silent=silent
;            endif
                extast,header,ast_info
                getrot,ast_info,rot,cdelt
                image_scale = 0.5*total(abs(cdelt)*3600.)
                ad2xy,ra_ref,dec_ref,ast_info,x_ref,y_ref
                image_size = size(image)

            ; check to see if there is any flux here
                kgphot,image,x_ref,y_ref,obj_flux,obj_sn,sky_flux,obj_flux_unc,sky_flux_unc, $
                  ap_rad=10.*image_scale,ap_sky=irac_ap_sky*image_scale,make_png=make_png, $
                  silent=silent,find_peak=0,image_scale=image_scale
                if (obj_flux EQ 0.0) then print,'***no good data ', files[i]
                if ((x_ref GE 0) AND (x_ref LT image_size[1]) AND $
                    (y_ref GE 0) AND (y_ref LT image_size[2]) AND $
                    (obj_flux NE 0.0)) then begin
                                ; bigger aperature to find the star
                    kgphot,image,x_ref,y_ref,obj_flux,obj_sn,sky_flux,obj_flux_unc,sky_flux_unc, $
                      ap_rad=10.*image_scale,ap_sky=irac_ap_sky*image_scale,make_png=make_png, $
                      silent=silent,find_peak=1,image_scale=image_scale
;                print,obj_flux
;                print,image_size
;                print,x_ref,y_ref

                    if (keyword_set(save_png)) then make_png = repstr(bfiles[z],'.fits','_phot') else make_png = 0
;            if (keyword_set(small_ap)) then begin ; find source
;                kgphot,save_image,x_ref,y_ref,obj_flux,obj_sn,sky_flux,obj_flux_un, $
;                       ap_rad=10.*image_scale,ap_sky=[10.,20.]*image_scale,make_png=make_png, $
;                       silent=silent,find_peak=find_peak,image_scale=image_scale
;            endif
                    if (keyword_set(test_centering)) then begin
                        save_x_ref0 = x_ref
                        save_y_ref0 = y_ref
                        kgphot,image,x_ref,y_ref,obj_flux,obj_sn,sky_flux,obj_flux_unc,sky_flux_unc, $
                          ap_rad=irac_ap_rad*image_scale,ap_sky=irac_ap_sky*image_scale,make_png=make_png, $
                          silent=silent,find_peak=find_peak,image_scale=image_scale
                        save_x_ref = x_ref
                        save_y_ref = y_ref
;                    print,x_ref,y_ref,obj_flux
                        
                        delta_x = 0.1
                        delta_y = 0.1
                        n_test = 21
                        max_flux = 0.0
                        for ii = 0,(n_test-1) do begin
                            x_ref = save_x_ref + delta_x*(ii - n_test/2)
                            for jj = 0,(n_test-1) do begin 
                                y_ref = save_y_ref + delta_y*(jj - n_test/2)
                                kgphot,image,x_ref,y_ref,obj_flux,obj_sn,sky_flux,obj_flux_unc,sky_flux_unc, $
                                  ap_rad=irac_ap_rad*image_scale,ap_sky=irac_ap_sky*image_scale,make_png=0, $
                                  silent=silent,find_peak=0,image_scale=image_scale
;                            print,obj_flux,format='(F8.4,$)'
                                if (max_flux LT obj_flux) then begin
                                    max_flux = obj_flux
                                    max_x_ref = x_ref
                                    max_y_ref = y_ref
                                endif
                            endfor
;                        print,' '
                        endfor
                        print,max_flux,max_x_ref,max_y_ref
                        print,max_x_ref-save_x_ref0,max_y_ref-save_y_ref0
                    
                                ; try with cntrd
                        x_ref = save_x_ref0
                        y_ref = save_y_ref0
                        cntrd,image,x_ref,y_ref,xcen,ycen,2.0
                        print,'***cntrd***'
                        print,x_ref,y_ref,xcen-x_ref,ycen-y_ref
                        print,xcen,ycen
                        print,'***cntrd***'
                        
;                    x_ref = max_x_ref
;                    y_ref = max_y_ref
                        x_ref = xcen
                        y_ref = ycen
                        find_peak = 0
                    endif else find_peak = 1
                    kgphot,image,x_ref,y_ref,obj_flux,obj_sn,sky_flux,obj_flux_unc,sky_flux_unc, $
                      ap_rad=irac_ap_rad*image_scale,ap_sky=irac_ap_sky*image_scale,make_png=make_png, $
                      silent=silent,find_peak=find_peak,image_scale=image_scale
                    ch_pos = strpos(files[i],'ch')
                    filter = 'IRAC'+strmid(files[i],ch_pos+2,1)
                    printf,unit1,names[k],filter,x_ref,y_ref,obj_flux,obj_flux/obj_sn,sky_flux,sky_flux_unc, $
                      format='(A15,A6,F8.2,F8.2,E12.4,E12.4,E12.4,E12.4)'
                    
                    for m = 0,(size_cube[3]-1) do begin
                        iimage = cube[*,*,m]
                        iimage_size = size(iimage)
                        iimage_scale = image_scale
                        ix_ref = x_ref
                        iy_ref = y_ref
                        if ((ix_ref GE 0) AND (ix_ref LT iimage_size[1]) AND $
                            (iy_ref GE 0) AND (iy_ref LT iimage_size[2])) then begin
                            if (keyword_set(save_png2)) then make_png = repstr(files[i],'.fits','_phot'+strtrim(string(m+1),2)) $
                            else make_png = 0
                            kgphot,iimage,ix_ref,iy_ref,obj_flux,obj_sn,sky_flux,obj_flux_un, $
                          ap_rad=irac_ap_rad*iimage_scale,ap_sky=irac_ap_sky*iimage_scale,make_png=make_png, $
                              silent=silent,find_peak=0,image_scale=iimage_scale
                            ch_pos = strpos(files[i],'ch')
                            filter = 'IRAC'+strmid(files[i],ch_pos+2,1)
                                ; correction factor for distortions
                            chn_num = fix(strmid(files[i],ch_pos+2,1)) - 1
;                            print,m, chn_num
                                ; correct the pixel positions given
                                ; these are subarrays
                            if (chn_num LT 2) then begin
                                ix_ref += 8
                                iy_ref += 216
                            endif else begin
                                ix_ref += 8
                                iy_ref += 8
                            endelse
                            corfac = a[0,chn_num] + a[1,chn_num]*(ix_ref - 128.) + a[2,chn_num]*(iy_ref - 128.) + $
                              a[3,chn_num]*(ix_ref - 128.)*(iy_ref - 128.) + a[4,chn_num]*(ix_ref - 128.)^2 + $
                              a[5,chn_num]*(iy_ref - 128.)^2
                                ; correction for pixel phase
                            if ((chn_num EQ 0) OR (chn_num EQ 1)) then begin
                                p_rad = sqrt(((ix_ref - fix(ix_ref)) - 0.5)^2 + ((iy_ref - fix(iy_ref)) - 0.5)^2)
                                corfac2 = 1.0 + irac_phase_a[chn_num]*((1./sqrt(2.*!PI)) - p_rad)
                            endif else corfac2 = 1.0
                                ; apply the correction
                            obj_flux *= corfac*corfac2
                            printf,unit2,names[k],m+1,filter,ix_ref,iy_ref,obj_flux,obj_flux/obj_sn,sky_flux,sky_flux_unc, $
                              chn_num,corfac,corfac2, $
                              format='(A15,I4,A6,F8.2,F8.2,E12.4,E12.4,E12.4,E12.4,I3,F8.2,F8.2)'
                        endif
                    endfor
                endif
            endfor
                                ; save the stacked image
            fits_write,repstr(files[i],'.lis','_med.fits'),output_image/n_bfiles,output_header

        endfor
    endif
endfor

free_lun,unit1
free_lun,unit2

end
