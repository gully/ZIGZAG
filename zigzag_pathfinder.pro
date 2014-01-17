;ZIGZAG: Zygo Interferometry of Gratings with Zernike Analysis by Gully

function read_zygo_xyz, xyz_fn, angle
;Author: gully
;Date: October 28, 2013
;Section I: function Read_Zygo_xyz
  ;Desc:     Reads in the raw .xyz file with no user intervention, saves as a .fits
  ;Inputs:   .xyz file, header info file
  ;Outputs:   .fits file with header info

;Outline:
;1) Read in the .xyz file header
;2) Read in the .xyz file data
;3) Turn the xyz data into a 2D array
;4) Add pars to the .fits header

;---------------------------------------------
;1) Read in the .xyz file header
readcol,xyz_fn, partnum, FORMAT = 'A', /SILENT, SKIPLINE = 6, numline=1
partnum_short=strcompress(partnum, /remove_all)
partnum_short=strmid(partnum_short, 1)
partnum_short=strcompress(partnum_short, /remove_all)
readcol,xyz_fn, src1, IntfScaleFactor1, WavelengthIn1, NumericAperture1, $
                ObliquityFactor1, Magnification1, CameraRes1, TimeStamp1, $
                FORMAT = 'I,F,F, I, I, I, F, D', /SILENT, SKIPLINE = 7, numline=1

wl_in=wavelengthIn1*1.0E6
cam_res=CameraRes1*1.0E6
caldat,linux2jd(timestamp1), Month, Day, Year, Hour, Minute, Second
time_string=strcompress(string(year)+string(month, format='(I02)')+string(day,format='(I02)'), /remove_all)
;---------------------------------------------


;---------------------------------------------
;2) Read in the .xyz file data
tag_arr= ['VERSION','DATASTART','DELIMITER','MISSINGVALUE','COMMENTSYMBOL',$
   'FIELDCOUNT','FIELDTYPES','FIELDNAMES','FIELDLOCATIONS','FIELDGROUPS']
tag_fmt='F,J,B,F,A,J(1),J(3),A(3),J(3),J(3)'
create_struct, at_new, '',tag_arr, tag_fmt  

at_new.version=1.0
at_new.datastart=14
at_new.delimiter=32
at_new.missingvalue=!Values.F_NaN
at_new.COMMENTSYMBOL=''
at_new.FIELDCOUNT=3
at_new.FIELDTYPES=[4,4,4]
at_new.FIELDNAMES=['FIELD1', 'FIELD2', 'FIELD3']
at_new.FIELDLOCATIONS=[0, 4, 7]
at_new.FIELDGROUPS=[0,0,0]

d0=read_ascii(xyz_fn, template=at_new)
x=reform(d0.field1[0, *])
y=reform(d0.field1[1, *])
z=reform(d0.field1[2, *])
n=n_elements(x)

x=x[0:n-2]
y=y[0:n-2]
z=z[0:n-2]
;---------------------------------------------


;---------------------------------------------
;3) Turn the xyz data into a 2D array
x1=max(x)*1.0
x0=min(x)*1.0
nx=x1-x0
y1=max(y)*1.0
y0=min(y)*1.0
ny=y1-y0

dx = (x1-x0)/nx
dy = (y1-y0)/ny
flat = fltarr(nx,ny)
flat[ (x-x0)/dx, (y-y0)/dy ] = z
;---------------------------------------------


;---------------------------------------------
;4) Add pars to the .fits header
mkhdr, hdr, flat
sxaddpar, hdr, 'CREATOR', 'gully', 'creator of this file'
sxaddpar, hdr, 'DATATYPE', 'Interferometry', 'Type of data- interferometry or PSF'
sxaddpar, hdr, 'TOOL', 'Zygo', 'Which tool was used'
sxaddpar, hdr, 'PARTNAME', partnum_short[0], 'name of the part'
sxaddpar, hdr, 'WAVELEN', wl_in[0], 'Wavelength of measurement in um'
sxaddpar, hdr, 'ANGLE', angle, 'angle of the measurement'
sxaddpar, hdr, 'PIXSCL', cam_res[0], 'Pixel scale (micron/pixel)'
sxaddpar, hdr, 'YMAX', ny*cam_res[0]/1000.0, 'y size of the beam in mm'
sxaddpar, hdr, 'XPROJ', nx*cam_res[0]/1000.0, 'Projected x size in mm'
sxaddpar, hdr, 'XDEPROJ', nx*cam_res[0]/1000.0/COS(angle/180.0*3.14159265), 'DeProjected x size in mm'
sxaddpar, hdr, 'ZUNITS', 'micron', 'units of z values'
sxaddpar, hdr, 'DATE', time_string[0], 'Date the data was taken, YYYYMMDD'
sxaddpar, hdr, 'ORIGFILE', xyz_fn, 'Original Filename'
sxaddpar, hdr, 'CLEANZZ', 0, 'ZigZag cleaning process (1 yes, 0 no)'
sxaddpar, hdr, 'ZSCALEZZ', 0, 'ZigZag zscale process (1 yes, 0 no)'
sxaddpar, hdr, 'ORIENTZZ', 0, 'ZigZag orient process (1 yes, 0 no)'
sxaddpar, hdr, 'DEPROJZZ', 0, 'ZigZag deproj process (1 yes, 0 no)'
sxaddpar, hdr, 'SPATFZZ', 0, 'ZigZag spatf process (1 yes, 0 no)'
sxaddpar, hdr, 'HPASS', 0, 'ZigZag spatf process (1 yes, 0 no)'
sxaddpar, hdr, 'LPASS', 0, 'ZigZag spatf process (1 yes, 0 no)'
sxaddpar, hdr, 'PASSKERN', 0.0, 'ZigZag spatf process (1 yes, 0 no)'
sxaddpar, hdr, 'SUBAPZZ', 0, 'ZigZag subap process (1 yes, 0 no)'
sxaddpar, hdr, 'APDMMZZ', 0, 'ZigZag subap diameter in mm'
sxaddpar, hdr, 'ZERNZZ', 0, 'ZigZag zern process (1 yes, 0 no)'
sxaddpar, hdr, 'ZERNTYZZ', 0, 'ZigZag zern type (FIT or RESIDUAL)'
sxaddpar, hdr, 'NZERNZZ', 0, 'ZigZag zern number of terms'
;sxaddhist,'Interferometry of G05 bottom part',HDR, /comment
;---------------------------------------------

outname=file_basename(xyz_fn, '.xyz')+'.fits'
writefits, outname, flat, HDR
return, 1
end

function read_zygo_fits, fits_name
;Aux1    function read_Zygo_fits
  ;Desc: reads the Zygo fits file and header, processes that info
  ;Inputs: fits file name
  ;Outputs: 2-D frame, header file in a structure
  
  d=readfits(fits_name, hdr)

  out_struct={HDR:hdr, $
              PARTNAME:sxpar(hdr, 'PARTNAME'), $
              WAVELEN:sxpar(hdr, 'WAVELEN'), $
              ANGLE:sxpar(hdr, 'ANGLE'), $
              PIXSCL:sxpar(hdr, 'PIXSCL'), $              
              YMAX:sxpar(hdr, 'YMAX'), $
              XPROJ:sxpar(hdr, 'XPROJ'), $
              XDEPROJ:sxpar(hdr, 'XDEPROJ'), $
              ZUNITS:sxpar(hdr, 'ZUNITS'), $
              CREATOR:sxpar(hdr, 'CREATOR'),$
              DATATYPE:sxpar(hdr, 'DATATYPE'),$
              TOOL:sxpar(hdr, 'TOOL'),$
              DATE:sxpar(hdr, 'DATE'),$
              ORIGFILE:sxpar(hdr, 'ORIGFILE'),$
              CLEANZZ:sxpar(hdr, 'CLEANZZ'),$
              ZSCALEZZ:sxpar(hdr, 'ZSCALEZZ'),$
              ORIENTZZ:sxpar(hdr, 'ORIENTZZ'),$
              SPATFZZ:sxpar(hdr, 'SPATFZZ'),$
              HPASS:sxpar(hdr, 'HPASS'),$
              LPASS:sxpar(hdr, 'LPASS'),$
              PASSKERN:sxpar(hdr, 'PASSKERN'),$    
              DEPROJZZ:sxpar(hdr, 'DEPROJZZ'),$
              SUBAPZZ:sxpar(hdr, 'SUBAPZZ'),$
              APDMMZZ:sxpar(hdr, 'APDMMZZ'),$
              ZERNZZ:sxpar(hdr, 'ZERNZZ'),$
              ZERNTYZZ:sxpar(hdr, 'ZERNTYZZ'),$
              NZERNZZ:sxpar(hdr, 'NZERNZZ')}

return, {dat:d, hd:out_struct}

end

function igram_clean, ds
;Section II: function igram_clean
  ;Desc:     Cleans the interferogram of bad pixels
  ;Inputs:   raw .fits file
  ;Outputs:  cleaned .fits
  
  swf=median(ds.dat, 10) 
  id=where(finite(ds.dat), complement=bi)
  ds.dat[bi]=swf[bi]
  
  hdr_out=ds.hd.hdr
  ds.hd.cleanzz=1
  sxaddpar, hdr_out, 'CLEANZZ', 1, 'ZigZag cleaning process (1 yes, 0 no)'
  ds.hd.hdr=hdr_out
  
  outname=file_basename(ds.hd.origfile, '.xyz')+'_ZZ_clean.fits'
  writefits, outname, ds.dat, hdr_out 
  
  return, ds

end

function igram_zscale, ds
;Section Ib: function igram_zscale
  ;Desc:     Scales the vertical values from phase to line edge position
  ;Inputs:   cleaned data
  ;Outputs:  z-scaled data with updated header
  
  ds.dat=ds.dat/sin(ds.hd.angle/180.0*3.141592654)

  hdr_out=ds.hd.hdr
  ds.hd.zscalezz=1
  sxaddpar, hdr_out, 'ZSCALEZZ', 1, 'ZigZag zscale process (1 yes, 0 no)'
  ds.hd.hdr=hdr_out
  
  outname=file_basename(ds.hd.origfile, '.xyz')+'_ZZ_zscale.fits'
  writefits, outname, ds.dat, hdr_out 


return, ds

end

function igram_orient, ds
;Section IIb: function igram_orient
  ;Desc:     Orient the interferogram-- rotate 180 and mirror.
  ;Inputs:   raw .fits file
  ;Outputs:  oriented .fits
  
  ds.dat=reverse(rotate(ds.dat, 2))

  hdr_out=ds.hd.hdr
  ds.hd.orientzz=1
  sxaddpar, hdr_out, 'ORIENTZZ', 1, 'ZigZag orient process (1 yes, 0 no)'
  ds.hd.hdr=hdr_out
  
  outname=file_basename(ds.hd.origfile, '.xyz')+'_ZZ_orient.fits'
  writefits, outname, ds.dat, hdr_out 


return, ds

end

function igram_deproj, ds
;Section III: function igram_deproj
  ;Desc:     Deprojects the interferogram from R3
  ;Inputs:   cleaned data
  ;Outputs:  deprojected .fits file with header info
  
  dat_dp=congrid(ds.dat, n_elements(ds.dat[*, 0])/cos(ds.hd.angle/180.0*3.141592654), n_elements(ds.dat[0, *]))
  
  hdr_out=ds.hd.hdr
  ds.hd.DEPROJZZ=1
  sxaddpar, hdr_out, 'DEPROJZZ', 1, 'ZigZag deproj process (1 yes, 0 no)'
  ds.hd.XPROJ=ds.hd.XDEPROJ
  sxaddpar, hdr_out, 'XPROJ', ds.hd.XDEPROJ
  ds.hd.hdr=hdr_out
  
  outname=file_basename(ds.hd.origfile, '.xyz')+'_ZZ_deproj.fits'
  writefits, outname, dat_dp, hdr_out 
  
  return, {dat:dat_dp, hd:ds.hd}

end  

function igram_spatf, ds
;Section IV: function igram_spatf
  ;Desc:     Spatial filter the input data
  ;Inputs:   data array
  ;Outputs:  smoothed data and residual data
  
  ;kernel
  kk=51 ; pixels (really this should be some fraction of the width
  smd=median(ds.dat, kk)
  res=ds.dat-smd
  
  hdr1=ds.hd.hdr
  hdr2=ds.hd.hdr
  sxaddpar, hdr1, 'SPATFZZ', 1, 'ZigZag spatf process (1 yes, 0 no)'
  sxaddpar, hdr2, 'SPATFZZ', 1, 'ZigZag spatf process (1 yes, 0 no)'
  sxaddpar, hdr1, 'HPASS', 1, 'ZigZag spatf process (1 yes, 0 no)'
  sxaddpar, hdr2, 'LPASS', 1, 'ZigZag spatf process (1 yes, 0 no)'
  sxaddpar, hdr1, 'PASSKERN', kk, 'ZigZag spatf process (1 yes, 0 no)'
  sxaddpar, hdr2, 'PASSKERN', kk, 'ZigZag spatf process (1 yes, 0 no)'

  
  out1=file_basename(ds.hd.origfile, '.xyz')+'_ZZ_HPASS.fits'
  out2=file_basename(ds.hd.origfile, '.xyz')+'_ZZ_LPASS.fits'
  writefits, out1, res, hdr1
  writefits, out2, smd, hdr2  
  
end  

function igram_subap, ds, cen, d_mm
;Section IV: function igram_subap
  ;Desc:     Crops a desired circular sub-aperture given a center and projected diameter
  ;Inputs:   data array, plate scake in micron per pixel, [xc, yc] array, diameter in mm
  ;Outputs:  cropped .fits file with header info
  
xc=cen[0]
yc=cen[1]

nx=n_elements(ds.dat[*, 0])
ny=n_elements(ds.dat[0, *])

d_pix=d_mm/ds.hd.pixscl*1000.0
r0=d_pix/2.0

if (xc+r0-1 gt nx) or (xc-r0 lt 0.0) or (yc+r0-1 gt ny) or (yc-r0 lt 0.0) then begin
  print, "-----------------------------------------------"
  print, "               ERROR                      "
  print, "Your desired beam diameter and/or center placement 
  print, "exceed the bounds of the observed interferometry"
  print, "  "
  print, "It looks like the beam size must be considerably less than:"
  print, ds.hd.pixscl/1000.0*nx < ds.hd.pixscl/1000.0*ny, 'mm'
  print, "-----------------------------------------------"
  stop
endif
orig_dat=ds.dat

r_arr=findgen(nx, ny)
mask_arr=r_arr*0.0+1.0
for i=0, nx-1 do begin
  for j=0, ny-1 do begin
    dx=(i-xc)
    dy=(j-yc)
    dr=sqrt(dx^2.0+dy^2.0)
    r_arr[i, j]=dr
    if dr gt r0 then mask_arr[i,j]=0
    if dr gt r0 then ds.dat[i,j]=65535 
  endfor
endfor

gx=where(ds.dat ne 65535) mod nx
gy=where(ds.dat ne 65535) / nx

sd=ds.dat[min(gx):max(gx), min(gy):max(gy)]


  hdr_out=ds.hd.hdr
  sxaddpar, hdr_out, 'YMAX', d_mm, 'y size of the beam in mm'
  sxaddpar, hdr_out, 'XPROJ', d_mm, 'Projected x size in mm'
  if ds.hd.DEPROJZZ eq 1 then $
      sxaddpar, hdr_out, 'XDEPROJ', d_mm, 'deProjected x size in mm'
  if ds.hd.DEPROJZZ eq 0 then $
      sxaddpar, hdr_out, 'XDEPROJ', d_mm/cos(ds.hd.angle/180.0*3.141592654), 'deProjected x size in mm'    
  sxaddpar, hdr_out, 'SUBAPZZ', 1, 'ZigZag subap process (1 yes, 0 no)'
  sxaddpar, hdr_out, 'APDMMZZ', d_mm, 'ZigZag subap diameter in mm'
  ds.hd.YMAX=d_mm
  ds.hd.XPROJ=d_mm
  ds.hd.XDEPROJ=d_mm/COS(ds.hd.angle/180.0*3.14159265)
  ds.hd.SUBAPZZ=1
  ds.hd.APDMMZZ=d_mm
  ds.hd.hdr=hdr_out
  
  outname=strcompress(file_basename(ds.hd.origfile, '.xyz')+'_'+$
          string(d_mm, format='(I2)')+'mm_ZZ_subap.fits', /remove_all)
  outname2=strcompress('ZZ_subap_mask_d'+string(d_mm, format='(I2)')+'mm.fits', /remove_all)
  writefits, outname, sd, hdr_out 
  mkhdr, hdr_mask, mask_arr
  sxaddpar, hdr_mask, 'D_MM', d_mm, 'Diameter of the ap. mask in mm'
  sxaddpar, hdr_mask, 'PIXSCL', ds.hd.PIXSCL, 'Pixel scale (micron/pixel)'
  writefits, outname2, (1.0-mask_arr)*orig_dat, hdr_mask

return, {dat:sd, hd:ds.hd}

end    

function igram_zern, ds, nz
;Section V: function igram_zern
  ;Desc:     Fits N zernike terms to the interferogram
  ;Inputs:   cleaned dat file, number of terms, partname, (flag to round up to filled radial shell?) 
  ;Outputs:  fitted .fits file, residual .fits file, coefficients in a .txt file.


device, decomposed=0
gif_set=1

wf0=ds.dat
; definitions
nx = n_elements(wf0[*, 0])    ; pixel in x-direction
ny = n_elements(wf0[0, *])    ; pixel in y-direction

; azimuthal frequency and radial order according to the Noll sequence
mn = NOLL_SEQUENCE( nz )

; Mask
mask = wf0
thresh=999.9
gi=where( mask lt thresh, complement=bi)
mask[bi]=0
mask[gi]=1
index = WHERE(mask, n_pix)
wf=wf0
wf[bi]=0.0

; Zernike polynomials as column vectors
zk = FLTARR( nz, n_pix )


FOR i = 1, nz - 1 DO BEGIN

    ; create a wavefront composed of Zernike polynomials
    tmp = FLOAT( ZERNIKE( mn[ 0, i ], mn[ 1, i ], nx ) )
    zk[ i, * ] = tmp[ index ]
    
ENDFOR


WINDOW, 10, XSIZE = 2 * nx, YSIZE = 2 * ny+70, XPOS=300, ypos=300
LOADCT, 13

; compute dynamic range for image display
m = MAX( ABS( wf[ index ] ) ) / 1.5 

; Singular Value Decomposition to get the coefficients
SVDC, zk, w, u, v
FOR i = 0, nz - 1 DO IF ABS( w[ i ] ) LE 1.0e-5 THEN w[ i ] = 0.0
cf = SVSOL( u, w, v, wf[ index ])

gd=0.0
FOR i = 0, nz - 1 DO BEGIN

    ; create a wavefront composed of Zernike polynomials
    tmp = FLOAT( ZERNIKE( mn[ 0, i ], mn[ 1, i ], nx ) )
    zk[ i, * ] = tmp[ index ]
    gd=gd+cf[i]* tmp * mask
    
    ; display the Zernike polynomials
    TV, BYTARR( 2*nx, 2 * ny+100 ), 0, 0
    this_zern=cf[i] * tmp * mask
    TV, BYTSCL( this_zern, min(this_zern), max(this_zern)), nx, ny
    TV, BYTSCL( gd * mask,  -1.0*m, 1.0*m ), 0, ny
    TV, BYTSCL( wf0* mask, -1.0*m, 1.0*m ), 0, 0
    resid=(wf0-gd)*mask
    TV, BYTSCL(resid, -3.0*stddev(resid), 3.0*stddev(resid)), nx, 0
    
    
    XYOUTS, 0, 2* ny+50, STRING( FORMAT = $
        '("Mode ordering number j: ", I3)', i ), $
        COLOR = 200, /DEVICE, FONT = 0
    XYOUTS, 0, 2* ny+30, STRING( FORMAT = $
        '("Azimuthal frequency m:  ", I3)', mn[ 0, i ] ), $
        COLOR = 200, /DEVICE, FONT = 0
    XYOUTS, 0, 2* ny+10, STRING( FORMAT = $
        '("Radial order n:         ", I3)', mn[ 1, i ] ), $
        COLOR = 200, /DEVICE, FONT = 0
    
    ;Labels     
    XYOUTS, 0, ny-10, ds.hd.partname+"!c measured!c wavefront", COLOR = 0, /DEVICE, FONT = 0
    XYOUTS, nx, ny-10, "Residual!c rescaled ", COLOR = 0, /DEVICE, FONT = 0
    XYOUTS, 0, 2*ny-10, "Net Zernike fit", COLOR = 0, /DEVICE, FONT = 0  
    XYOUTS, nx, 2*ny-10, "jth Zernike !c rescaled ", COLOR = 0, /DEVICE, FONT = 0
    
    WAIT, 0.1

  if gif_set then begin    
      WSET, 10
      image24=TVRD(true=1)
      image2d=color_quan(image24, 1, r, g, b)
      outfile=strcompress(ds.hd.partname+"_tempgifzern_"+string(i,FORMAT="(I03)"), /remove_all)
      write_gif, outfile+'.gif', image2d, r, g, b
  endif
   
    
ENDFOR

spawn, 'gifsicle --delay=10 --loop *_tempgifzern_* > '+ $
        strcompress(ds.hd.partname+'_Zern_anim.gif', /remove_all)
spawn, 'rm *_tempgifzern_*'

!P.MULTI = 0

residual=wf0-gd
residual[bi]=65535
gd[bi]=65535

hdr1=ds.hd.hdr
hdr2=ds.hd.hdr
  
sxaddpar, hdr1, 'ZERNZZ', 1, 'ZigZag zern process (1 yes, 0 no)'
sxaddpar, hdr2, 'ZERNZZ', 1, 'ZigZag zern process (1 yes, 0 no)'
sxaddpar, hdr1, 'ZERNTYZZ', 'FIT', 'ZigZag zern type (FIT or RESIDUAL)'
sxaddpar, hdr2, 'ZERNTYZZ', 'RESIDUAL', 'ZigZag zern type (FIT or RESIDUAL)'
sxaddpar, hdr1, 'NZERNZZ', nz, 'ZigZag zern number of terms'
sxaddpar, hdr2, 'NZERNZZ', nz, 'ZigZag zern number of terms'
  
out1=file_basename(ds.hd.origfile, '.xyz')+'_ZZ_fit.fits'
out2=file_basename(ds.hd.origfile, '.xyz')+'_ZZ_resid.fits'
out3=file_basename(ds.hd.origfile, '.xyz')+'_ZZ_Zcoeffs.txt'
com3=strcompress(';index, Az #, Rad #, coeff ('+ds.hd.ZUNITS+')')
writefits, out1, gd, hdr1
writefits, out2, residual, hdr2 
forprint, indgen(nz), mn[ 0, * ], mn[ 1, * ], cf, textout=out3, COMMENT=com3

return, 1

end  

function igram_psf, ds, par
;Section X: function igram_orient
  ;Desc:     make a PSF given a circularized sub-aperture
  ;Inputs:   sup-apertured .fits file
  ;Outputs:  psf file.
  
  d=ds.dat
  bi=where(d eq 65535, complement=gi)
  
  
  d[bi]=par
  nx=n_elements(d[*, 0])
  ny=n_elements(d[0, *])
  ef=6.0
  dex=fltarr(nx*ef, nx*ef)+par
  dex[nx*ef/2.0:nx*ef/2.0+nx-1, nx*ef/2.0:nx*ef/2.0+ny-1]=d
  ftd=fft(dex)
  sfd=shift(ftd, nx*ef/2.0-1.0, nx*ef/2.0-1.0)
  psd=abs(sfd)^2.0
  
  psd=psd/max(psd)
  psc=congrid(psd, nx, ny)
  psc=psc/max(psc)
  ;tvscl, bytscl(alog10(psd), -6, 0)
  writefits, 'TEMP_PSF.fits', psd
  
  ;Zygo data
  alt_dir='/Volumes/cambridge/Astronomy/silicon/CA1/interferom/20120604_gully/gully/'
  restore, filename=alt_dir+'zygo_PSF_file_ascitemplate.sav'
  d2=read_ascii('E04_R3_25mm.psf', template=at)
  dat2=d2.field001
  dat2=dat2*dat2
  dat2=dat2/max(dat2)
  
  n_lam_d=20.0
  res_el=3.0
  mxval=max(psd, id)
  idx=id mod fix(nx*ef)
  idy=id / fix(nx*ef)
  
  window, 10, xsize=256, ysize=256, xpos=400, ypos=0
  window, 11, xsize=256, ysize=256, xpos=700, ypos=0
  window, 12, xsize=256, ysize=256, xpos=100, ypos=0
  
  
;  subp=psd[idx-res_el*n_lam_d:idx+res_el*n_lam_d, $
;            idy-res_el*n_lam_d:idx+res_el*n_lam_d]

  subp=psd[idx-128:idx+128, idy-128:idy+128]
  wset, 11
  tvscl, bytscl(alog10(subp), -6, 0)
  wset, 10
  tvscl, bytscl(alog10(dat2), -6, 0)
  wset, 12
  tvscl, bytscl(d, -1,1)
  
  print, 1
  
  return, 1
  
end

pro zigzag_pathfinder
;Author: gully
;date: October 28, 2013
;desc: 
;
;ZIGZAG: Zygo Interferometry of Gratings with Zernike Analysis by Gully

;Section 0: inputs
device, RETAIN=2
cd, '/Volumes/cambridge/Astronomy/silicon/Heritage/E04/Zygo/'
xyz_fn='E04_R3.xyz'
fits_name='E04_R3.fits'  
angle=71.6 ;degrees
d_mm=25.0 ;mm
cen=[155, 400] ;pixels
nz=91


;Section I: function Read_Zygo_xyz
out1=read_zygo_xyz(xyz_fn, angle)

ds=read_zygo_fits(fits_name)

;Section II: function igram_clean
d_cln=igram_clean(ds)

;--------------------------------------------
;Below here we make ENGINEERING data products
;--------------------------------------------

;Section IIb: function igram_clean
d_zs=igram_zscale(d_cln)

;Section IIc: function igram_clean
d_ori=igram_orient(d_zs)

;Section III: function igram_deproj
d_dp=igram_deproj(d_ori)

;Section IVa: function igram_filt
out1=igram_spatf(d_dp)


;--------------------------------------------
;From here on we make ASTRONOMY data products
;--------------------------------------------

;Section IV: function igram_subap
d_sa=igram_subap(d_cln, cen, d_mm)

;Section V: function igram_zern
out2=igram_zern(d_sa, nz)

;Section V: function igram_PSF

d_sa.dat=readfits('E04_R3_ZZ_resid.fits')

;d_sa.dat=readfits('E04_R3_ZZ_fit.fits')
out3=igram_psf(d_sa, 0.0)


;TODO: PSF

print, 1

end