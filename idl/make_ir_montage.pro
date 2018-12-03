; program to make a montage of the IRAC, IRS-blue, and MIPS 24um data
; for a single star

pro make_ir_montage,name

tname = repstr(name,'d0','d')

irac_files = file_search('~/Spitzer/IR_Ext/IRAC/data/' + tname + '*.fits',count=n_irac_files)

irs_files = file_search('~/Spitzer/IR_Ext/IRS_peakup/data/' + tname + '*.fits',count=n_irs_files)

readcol,'~/Dust/Ext/MIPS/irext_mips24_all.dat', $
        format='(A5,L12,A15,F15.3,F12.3,F10.3,F12.3,F12.3,F10.3,F12.3)', $
        tag,aorkey,mname,obstime,aflux,asn,afluxmjy,pflux,psn,flux,/silent

n_names = n_elements(mname)
newname = strarr(n_names)
for i = 0,(n_names-1) do begin
    newname[i] = strlowcase(mname[i])
    if (strlen(newname[i]) EQ 7) then newname[i] = repstr(newname[i],'hd','hd0')
    if (newname[i] EQ 'ngc2024') then newname[i] += '_1'
    if (newname[i] EQ 'bd+631964') then newname[i] = 'bd+63d1964'
endfor

mindxs = where(name EQ newname,n_mindxs)
if (n_mindxs GT 0) then begin
    mips_files = file_search('~/Spitzer/IR_Ext/MIPS/data/*' + strtrim(aorkey[mindxs],2) + '*_mos.fits',count=n_mips_files)

    if (n_mips_files EQ 0) then begin
        mips_files = file_search('~/Spitzer/IR_Ext/MIPS/data/*' + strtrim(strlowcase(mname[mindxs]),2) + '*_mos.fits',count=n_mips_files)
    endif
endif else begin
    mips_files = ''
    n_mips_files = 0
endelse

print,'IRAC = ', irac_files
print,'IRS = ', irs_files
print,'MIPS = ', mips_files

; read the files

if (n_mips_files GT 0) then begin
    fits_read,mips_files[0],mips_image,mips_header 
    getrot,mips_header,mips_rot,mips_cdelt
    mips_pixscale = 0.5*total(abs(mips_cdelt))*3600.
endif else begin 
    mips_image = fltarr(182,199)
    mips_pixscale = 2.49
endelse
if (n_irs_files GT 0) then begin
    fits_read,irs_files[0],irs_image,irs_header
    getrot,irs_header,irs_rot,irs_cdelt
    irs_pixscale = 0.5*total(abs(irs_cdelt))*3600.
endif else begin
    irs_image = fltarr(84,114)
    irs_pixscale = 1.83255
endelse

fits_read,irac_files[0],irac_image1,irac_header1
fits_read,irac_files[1],irac_image2,irac_header2
fits_read,irac_files[2],irac_image3,irac_header3
fits_read,irac_files[3],irac_image4,irac_header4

; determine the size of the images
getrot,irac_header1,irac_rot,irac1_cdelt
getrot,irac_header2,irac_rot,irac2_cdelt
getrot,irac_header3,irac_rot,irac3_cdelt
getrot,irac_header4,irac_rot,irac4_cdelt

irac1_pixscale = 0.5*total(abs(irac1_cdelt))*3600.
irac2_pixscale = 0.5*total(abs(irac2_cdelt))*3600.
irac3_pixscale = 0.5*total(abs(irac3_cdelt))*3600.
irac4_pixscale = 0.5*total(abs(irac4_cdelt))*3600.

print,'MIPS pscale = ', mips_pixscale
print,'IRS pscale = ', irs_pixscale
print,'IRAC1 pscale = ', irac1_pixscale
print,'IRAC2 pscale = ', irac1_pixscale
print,'IRAC3 pscale = ', irac1_pixscale
print,'IRAC4 pscale = ', irac1_pixscale

mips_size = fix(round(mips_pixscale*(size(mips_image))[1:2]))
irs_size = fix(round(irs_pixscale*(size(irs_image))[1:2]))
irac1_size = fix(round(irac1_pixscale*(size(irac_image1))[1:2]))
irac2_size = fix(round(irac2_pixscale*(size(irac_image2))[1:2]))
irac3_size = fix(round(irac3_pixscale*(size(irac_image3))[1:2]))
irac4_size = fix(round(irac4_pixscale*(size(irac_image4))[1:2]))

print,'MIPS size = ', mips_size
print,'IRS size = ', irs_size
print,'IRAC1 size = ', irac1_size
print,'IRAC2 size = ', irac2_size
print,'IRAC3 size = ', irac3_size
print,'IRAC4 size = ', irac4_size

; scale the individual images to run from 0 to 1

sub_mips_image = mips_image[75:125,75:125]
indxs = where(finite(sub_mips_image))
std_mips = stdev(sub_mips_image[indxs],mean_mips)
mips_image = (mips_image - mean_mips + std_mips)/(3.0*std_mips)

indxs = where(finite(irs_image))
std_irs = stdev(irs_image[indxs],mean_irs)
irs_image = (irs_image - mean_irs + std_irs)/(3.0*std_irs)

indxs = where(finite(irac_image1))
std_irac = stdev(irac_image1[indxs],mean_irac)
irac_image1 = (irac_image1 - mean_irac + std_irac)/(3.0*std_irac)

indxs = where(finite(irac_image2))
std_irac = stdev(irac_image2[indxs],mean_irac)
irac_image2 = (irac_image2 - mean_irac + std_irac)/(3.0*std_irac)

indxs = where(finite(irac_image3))
std_irac = stdev(irac_image3[indxs],mean_irac)
irac_image3 = (irac_image3 - mean_irac + std_irac)/(3.0*std_irac)

indxs = where(finite(irac_image4))
std_irac = stdev(irac_image4[indxs],mean_irac)
irac_image4 = (irac_image4 - mean_irac + std_irac)/(3.0*std_irac)

; resample everything to 1"/pixel

mips_image = congrid(mips_image,mips_size[0],mips_size[1])
irs_image = congrid(irs_image,irs_size[0],irs_size[1])
irac_image1 = congrid(irac_image1,irac1_size[0],irac1_size[1])
irac_image2 = congrid(irac_image2,irac2_size[0],irac2_size[1])
irac_image3 = congrid(irac_image3,irac3_size[0],irac3_size[1])
irac_image4 = congrid(irac_image4,irac4_size[0],irac4_size[1])

; new image

nimage_size = [mips_size[0]+irs_size[0],mips_size[1]]
nimage = replicate(!values.f_nan,nimage_size[0],nimage_size[1])

; fill the image

nimage[irs_size[0]:nimage_size[0]-1,*] = mips_image
nimage[0:irs_size[0]-1,40:irs_size[1]-1+40] = irs_image
nimage[20:irac1_size[0]-1+20,nimage_size[1]-irac1_size[1]-20:nimage_size[1]-1-20] = irac_image1
nimage[20:irac2_size[0]-1+20,nimage_size[1]-irac1_size[1]-irac2_size[1]-40:nimage_size[1]-irac1_size[1]-1-40] = irac_image2
nimage[20:irac3_size[0]-1+20,nimage_size[1]-irac1_size[1]-irac2_size[1]-irac3_size[1]-60:nimage_size[1]-irac1_size[1]-irac2_size[1]-1-60] = irac_image3
nimage[20:irac4_size[0]-1+20,nimage_size[1]-irac1_size[1]-irac2_size[1]-irac3_size[1]-irac4_size[1]-80:nimage_size[1]-irac1_size[1]-irac2_size[1]-irac3_size[1]-1-80] = irac_image4

; save the image

fitsname = 'ir_ext/' + name + '_montage.fits'
fits_write,fitsname,nimage

; save a png version

create_png,fitsname,repstr(fitsname,'.fits','.png'),int_range=[0.,1.5],color_table=0,/reverse

end
