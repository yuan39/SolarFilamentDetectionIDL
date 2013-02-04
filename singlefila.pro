Pro singleFila,filepath,outpath
;outpath='D:\DisappearanceList\'
;filepath = 'D:\DisappearanceList\kanz_halph_fd_20100801_0847.hi.fts.gz'

loadct,0

img_ori = readfits(filepath,header)

xcen = round(sxpar(header,'CENTER_X'))
ycen = round(sxpar(header,'CENTER_Y'))
rad = round(sxpar(header,'SOLAR_R'))
datestr = sxpar(header,'DATE')
solar_p0 = round(sxpar(header,'SOLAR_P0'))

img = img_ori(xcen-rad:xcen+rad,ycen-rad:ycen+rad)
img = rot(img,solar_p0)

clean_pixel_number = 200

sz = size(img,/dimension)
width = sz[0]
height = sz[1]

mask = findgen(width,height)*0.0
mask[where(shift(dist(width,height),fix(width/2),fix(height/2)) le width/2*0.9) ] = 1;

mean_img = mean(img[where(mask gt 0)])
std_img = stdev(img[where(mask gt 0)])

print,mean_img
print,std_img

steps = 13
ts = findgen(steps)
for i = 0,steps-1 do begin
  ts[i] = mean_img-3*std_img + abs(3*std_img)/steps*(i+1)
endfor

kernel = replicate(1,3,3)
kernel[1,1]=-8

img_non_mask = img;
img = img * mask

g = convol(img,kernel,center=1)

thresh = mean_img-3*std_img
tmp_tmp = 0.0
;forcurve = findgen(steps-1)
for i = 0,steps-2 do begin
  img1 = img lt ts[i]
  img2 = img lt ts[i+1]
  diff_region = img2 gt img1;
  if total(diff_region) gt 1.0 then begin
    if total(diff_region*g)/total(diff_region) gt tmp_tmp then begin
      thresh = ts[i]
      tmp_tmp = total(diff_region*g)/total(diff_region)
    endif
  endif
  
  ;forcurve(i) = total(diff_region*g)/total(diff_region)
  ;tvscl,congrid((img lt ts [1]) gt (img lt ts[0]),600,600)
  ;window,i+1,xs=600,ys=600
  ;tvscl,congrid(diff_region,600,600)
endfor  
;window,10
;plot,ts(0:steps-1),forcurve

;thresh =  ts[round(mean(where(forcurve eq max(forcurve))))]
threshed = (img lt thresh)*mask;
window,11,xs=600,ys=600
tvscl,congrid(threshed,600,600)

labeled = label_region(threshed,/ALL_NEIGHBORS, /ULONG)
hist = histogram(labeled)
cleaned = findgen(width,height)*0.0
for i=1,n_elements(hist)-1 do begin
  if hist[i] gt clean_pixel_number then begin
    cleaned[where(labeled eq i)] = 1.0
  endif
endfor
cleaned = cleaned gt 0.0
window,12,xs=600,ys=600
tvscl,congrid(cleaned,600,600)

se_size = 50
se = indgen(se_size,se_size)*0
se[where(shift(dist(se_size,se_size),fix(se_size/2),fix(se_size/2)) le se_size/2) ] = 1;
closing = MORPH_CLOSE(cleaned, se)
window,13,xs=600,ys=600
tvscl,congrid(closing,600,600)




;;;;;;;;;;;;;;;;;final clean%%%%%%%%%%%%%%%%%%

closing_labeled = label_region(closing,/ALL_NEIGHBORS, /ULONG)
hist = histogram(closing_labeled)
final = findgen(width,height)*0.0
for i=1,n_elements(hist)-1 do begin
  if hist[i] gt 100 then begin
    final[where(closing_labeled eq i)] = 1.0
  endif
endfor
final = final gt 0.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;remove sunspots;;;;;;;;;;;;;



dims = SIZE(final, /DIMENSIONS) 
print,dims
; Get blob indices:
final_b = LABEL_REGION(final)
; Get population and members of each blob:
h = HISTOGRAM(final_b, REVERSE_INDICES=r)
FOR i=1, N_ELEMENTS(h)-1 DO BEGIN
   ;Find subscripts of members of region i.
   p = r[r[i]:r[i+1]-1]  
   
   perimeterIndex = Find_Boundary(p,  XSize=dims[0], YSize=dims[1],Perimeter=length, Area=area)
          
   Print, 'Length: ', length, '     Area: ', area
   metric = 4*3.1415926*area/length^2
   print,metric
   if (metric ge 0.5 ) THEN BEGIN
      final[p]=0
   endif
   
ENDFOR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
thisDevice = !D.Name
Set_Plot,'Z'
Erase

dims = SIZE(final, /DIMENSIONS) 
Device, Set_Resolution=[dims[0],dims[1]],Set_Pixel_Depth = 24, Decomposed=1
tvimage,bytscl(img_non_mask)

final_labeled = label_region(final,/ALL_NEIGHBORS, /ULONG)
dims = SIZE(final_labeled, /DIMENSIONS) 
hist = HISTOGRAM(final_labeled, REVERSE_INDICES=r)

lats = findgen(n_elements(hist)-1)
lons = findgen(n_elements(hist)-1)
areas = findgen(n_elements(hist)-1)

openw, 1, outpath+datestr+'.result.txt'
printf,1,'% '+ datestr
print,'% '+datestr
printf,1,'% id,Area (square megameter),longitude (degrees),latitude (degrees)'
for i=1,n_elements(hist)-1 do begin
   
  index = where(final_labeled eq i)
  xcoor = 0;
  ycoor = 0;
  for j=0,n_elements(index)-1 do begin
      xcoor = xcoor + index[j] mod width
      ycoor = ycoor + index[j] / width
  endfor
  xcoor = xcoor/hist[i]
  ycoor = ycoor/hist[i]
  
  lats[i-1]= asin(float(ycoor - float(height/2))/rad)*180/!PI
  tmp_r = sqrt(rad*rad - (ycoor - float(height/2))^2)
  lons[i-1] = asin( float(xcoor - float(width/2))/tmp_r)*180/!PI
  areas[i-1] = hist[i] * (1392000.0/width/1000.0)^2 
 
  ;PLOTS, CIRCLE(fix(float(xcoor)/width*600), fix(float(ycoor)/height*600), 10 ),Color=FSC_Color('Red'), /Device
  
  p = r[r[i]:r[i+1]-1]  
  PLOTS, Find_Boundary(p, XSize=dims[0], YSize=dims[1], Perimeter=length, Area=area), $
           /Device, Color=FSC_Color('red'),thick=2

  ;XYOuts, fix(float(xcoor)/width*600)+10, fix(float(ycoor)/height*600), /Device, STRING(i, FORMAT='(I2.2)'), Color=FSC_Color('Red')
  XYOuts, fix(float(xcoor))+10, fix(float(ycoor)),charsize=6, /Device, STRING(i, FORMAT='(I2.2)'), Color=FSC_Color('Red'),charthick=4
  
  
  
  printf, 1, i,areas[i-1],lons[i-1],lats[i-1],format='(i,",",f0.1,",",f0.1,",",f0.1)'
  print,  i,areas[i-1],lons[i-1],lats[i-1],format='(i,",",f0.1,",",f0.1,",",f0.1)'
  
endfor
close,1

xyouts,60,height-80,"N",charsize=6,charthick=4,/Device
xyouts,0,height-160,"E  W",charsize=6,charthick=4,/Device
xyouts,60,height-240,"S",charsize=6,charthick=4,/Device
xyouts, 10, 10, "KANZ",charsize=6,charthick=4,/Device
xyouts,width/2-80,10,sxpar(header,'DATE_OBS') , charsize=6,charthick=4,/Device


snapshot = TVRD(True=1)
Set_Plot, thisDevice
;write_jpeg, outpath+datestr+'.char.jpg',snapshot, True=1, Quality = 90 

write_png,outpath+datestr+'.char.png', snapshot 


print,outpath+datestr+'.result.txt'
print,outpath+datestr+'.char.png'
file_copy,outpath+datestr+'.result.txt','/var/www/current/result.txt',/force,/overwrite
file_copy,outpath+datestr+'.char.png','/var/www/current/char.png',/force,/overwrite
end


