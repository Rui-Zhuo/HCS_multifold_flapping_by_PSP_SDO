read_sdo,file,hdr,image
image=image/hdr.exptime;normalize exposure time
imdisp=mgn(image,h=0.95)
tv,hist_equal(imdisp,percent=0.05)

imdisp=mgn(image,h=0.925)
indnan=where(~finite(imdisp),cntnan,comp=nindnan);identify NANs
imdisp=(imdisp>minval)<maxval      ;threshold output
imdisp=(imdisp-minval)*255/(maxval-minval)    ; normalize to 255 max
imdisp=byte(imdisp)     ;convert to byte
if cntnan gt 0 then imdisp[indnan]=200      ;set NANs to 200 (light grey)

End