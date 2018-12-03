
pro get_mips24_from_scan,filename,name,coords

fits_open,filename,fcb
fits_open,'data/'+name+'_mips24_from_scan_mos.fits',ofcb,/write

for i = 0,3 do begin
    if (i EQ 0) then pdu = 1 else pdu = 0
    fits_read,fcb,image,header,pdu=pdu,exten_no=(i+1)

    new_header = header
    sxaddpar,new_header,'NAXIS1',182
    sxaddpar,new_header,'NAXIS2',199
    sxaddpar,new_header,'CRVAL1',coords[0]
    sxaddpar,new_header,'CRVAL2',coords[1]
    sxaddpar,new_header,'CRPIX1',91.5
    sxaddpar,new_header,'CRPIX2',97.5
    sxaddpar,new_header,'RA_REF',coords[0]
    sxaddpar,new_header,'DEC_REF',coords[1]

    hastrom,image,header,new_header
    
    fits_write,ofcb,image,new_header
endfor

end
