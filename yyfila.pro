Pro yyfila,daystr = daystr, monthstr = monthstr, yearstr = yearstr

if not keyword_set(daystr) then begin

  caldat,systime(/julian),month,day,year
  daystr = STRING(day, FORMAT='(I2.2)')
  monthstr = STRING(month, FORMAT='(I2.2)')
  yearstr = STRING(year, FORMAT='(I4.4)')

endif

files = sock_find('http://cesar.kso.ac.at','kanz_halph_fc_'+yearstr+monthstr+daystr+'*fts.gz',path='/halpha4M/FITS/normal/'+yearstr)

if (size(files,/dimension) gt 0) then begin
  filepath = '/var/www/halpha/'+yearstr+monthstr+daystr+'.fts.gz'
  sock_copy,files[0],filepath,/verb,/outdir
endif else begin
  print,'no file for '+ yearstr+monthstr+daystr
  return
endelse

singleFila,filepath,'/var/www/current/'

end
