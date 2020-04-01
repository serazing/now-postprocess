;Name:
;       calculate_mslp_arguments.ncl
;
;Purpose:
;       This script reads various variables (3-d T, surface pressure, geopotential height)
;       from wrfout files and calculates the mean sea level pressure field using
;       the ncl wrf_user_getvar(f[ff],"slp",-1) function.
;       It then writes the mslp field in a new netcdf file.

;Calling sequence:
;       ./runncl.sh calculate_mslp_arguments.ncl arg1 arg2 arg3 arg4 
;
;Author:
;       Alejandro DI Luca
;       a.diluca@unsw.edu.au
;       Based on Fei Ji's script (slp_3hrly.ncl)
;
;Creation:
;       05/07/2013
;
;Comments:
;     18-11-2013:
;     - I modified output directory form and added several variables to specify 
;       the version of the RCM and the domain.  
;     8-07-2014:
;     - I modified added two new arguments, year and month, so now the python script 
;       used to give the arguments can check if output files already exists.
;
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

; set the number arguments

if(ismissing(getenv("NCL_N_ARGS"))) then 
   NCL_N_ARGS = 0  
else
   NCL_N_ARGS = getenv("NCL_N_ARGS")
end if

; test for correct number of arguments (optional)
if (NCL_N_ARGS .lt. 4) then
  print("This script needs at least 2 command line arguments. Exit.")
  print("Usage: ./runncl.sh runwithargs.ncl $YEAR")
  exit()
end if

; read arguments
name_script = getenv("NCL_ARG_1")
print("First Argument (name script) ="+name_script)

indir = getenv("NCL_ARG_2")
print("Second Argument (input dir) ="+indir)

outdir = getenv("NCL_ARG_3")
print("Third Argument (output dir) ="+outdir)

year = getenv("NCL_ARG_4")
print("4th Argument (year) ="+year)

month = getenv("NCL_ARG_5")
print("5th Argument (month) ="+month)

strsplit=str_split(outdir,"/")
print((/strsplit/))

; Get gcm, rcm, year, domain names
domain_name="d01"

file_head="wrfout_"+domain_name+"_"

;get the 2 digit month string
monstr = sprinti("%0.2i",stringtointeger(month))
     
; get the files of interest
files = systemfunc("ls "+indir+file_head+year+"-"+monstr+"-*")
print((/files(0)/))
     
;load the files
f = addfiles(files+".nc","r")
     
; Read Times variable
Times = f[:]->Times
     
frec = dimsizes(Times(:,0))

;the last file each month has the 0 hour timestep for the next month
;the first file each month also has the 0 hour
;check that this is the case and only count one
  
;check the 0 hour is the last time step in the file
;if it is don't include it
if (isStrSubset(chartostring(Times(frec-1,:)),"-01_00:00:00")) then
      zeroHour = True
      frec = frec - 1
else
      zeroHour = False
end if

newTimes = Times(:frec-1,:)     
    
; Get the lat & lon arrays if I haven't already
if (.not.isvar("lat2d")) then
lat2d = f[0]->PSFC
end if
    
    ; Creates output SLP variable for this month
    slp = new((/frec,dimsizes(lat2d(0,:,0)),dimsizes(lat2d(0,0,:))/),"float")
    copy_VarCoords(newTimes(:,0),slp(:,0,0))
    copy_VarCoords(lat2d,slp(0,:,:))
    slp@MemoryOrder = "XY"
    slp@FieldType = "104"
    slp@coordinates = "XLONG XLAT"
    slp@history = "ncl wrf_user_getvar(slp)"
    slp@description = "Sea Level Pressure"
    slp@units = "hPa"

;   ######################################################
    ; FILES LOOP
    rec = 0	
    do ff = 0,dimsizes(files)-1
      print(ff(0))

      ; Get slp from WRF output
      var = wrf_user_getvar(f[ff],"slp",-1)
      
      if (ff.eq.dimsizes(files)-1.AND.zeroHour) then
        erec = dimsizes(var(:,0,0))-2
        slp(rec:,:,:) = (/var(:erec,:,:)/)
      else
        erec = dimsizes(var(:,0,0))-1
        slp(rec:rec+erec,:,:) = (/var(:,:,:)/)
      end if
      
      rec = rec + erec + 1
      
      delete(var)
      
    end do	; file's loop
    
    ;write output
    system("rm "+outdir+"WRF_mslp_"+domain_name+"_"+year+"-"+monstr+".nc")
    fout = addfile(outdir+"WRF_mslp_"+domain_name+"_"+year+"-"+monstr+".nc","c")
    print (outdir+"WRF_mslp_"+domain_name+"_"+year+"-"+monstr+".nc")

    ; Make time an UNLIMITED dimension
    filedimdef(fout,"Time",-1,True)

    fout->SLP = slp
    
    fout->Times = newTimes
    
    ; Delete variables
    delete(slp)
    delete(Times)
    delete(newTimes)
    delete(files)
    delete(f)
 
end
