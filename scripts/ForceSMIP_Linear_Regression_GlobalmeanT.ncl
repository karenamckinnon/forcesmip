; Estimates the forced response of a single run by regressing on global-mean temperature. 
; 
; This script calls a python script (python butterworthfilter.aa.py) to apply a 10yr lowpass 
; butterworth filter. Before running this script on NCAR's casper machine, enter the following
; in the terminal so a proper python environment (and NCL) is loaded:
;   module load ncarenv conda
;   conda activate npl-2023a
;   module load ncl
; If upon running this script you get error messages that state: 
;   ImportError: No module named xarray, warning:_NclOpenFile: cannot open file <tmp2.nc>; No such file or directory
; then your python environment is not set up correctly. Try entering: 
;   module purge 
; and then again enter the the three lines above.
;
; The only user inputs are:
; - The file name (line 30) 
; - The start year of the analysis (line 31)
; - The end year of the analysis (line 32)
;
; This script may be run on any file in the Training, Evaluation-Tier1, Evaluation-Tier2, or Evaluation-Tier3 
; directories. The only exeption is that you may not specify a file that contains the zmta variable.
; 
; To run this script, set lines 30-32 as desired, then in the terminal enter:
;   ncl ForceSMIP_Linear_Regression_GlobalmeanT.ncl 
;
begin
;  filein is the file to be analyzed
;    Example Evaluation file setting: filein = "/glade/campaign/cgd/cas/asphilli/ForceSMIP/Evaluation-Tier3/Amon/tas/tas_mon_3I.197901-202212.nc"
;    Example Training file setting: filein = "/glade/campaign/cgd/cas/asphilli/ForceSMIP/Training/Aday/monmaxtasmax/MPI-ESM1-2-LR/monmaxtasmax_day_MPI-ESM1-2-LR_historical_ssp585_r9i1p1f1.188001-202212.nc"
  filein = "/glade/campaign/cgd/cas/asphilli/ForceSMIP/Evaluation-Tier1/Amon/pr/pr_mon_1I.195001-202212.nc"
  syear = 1955    ; start year of analysis
  eyear = 2022    ; end year of analysis
;-------------------------------------
;  Parse the file name for variable name, start date and end date of files
;
  temp := str_split(filein,"/")
  fn = temp(dimsizes(temp)-1)

  tm = str_split(fn,".")
  ttt = str_split(tm(dimsizes(tm)-2),"-")
  syearF = toint(ttt(0))   ; start year/start month
  eyearF = toint(ttt(1))   ; end year/end month

  vT = str_split(fn,"_")
  varT = vT(0)
;-------------------------------------
; Figure out the name of the matching tas file by parsing the given file.
;
  filein_tas = "/"+str_join(temp(:6),"/")    
  if (temp(6).eq."Training") then
     filein_tas = filein_tas+"/Amon/tas/"+temp(9)+"/"+str_sub_str(fn,varT,"tas")
  else
     filein_tas = filein_tas+"/Amon/tas/"+str_sub_str(fn,varT,"tas")
  end if
  filein_tas = str_sub_str(filein_tas,"tas_day_","tas_mon_")
  print("tas file identified as "+filein_tas)
;-------------------------------------
;  Check to make sure that zmta is not the variable to be analyzed, and that the specified years to be analyzation fall within syearF/eyearF
;
  if (varT.eq."zmta") then
     print("error: This script may not be used to analyze the zmta variable, exiting")
     exit
  end if
  if (syear.lt.(syearF/100).or.eyear.gt.(eyearF/100)) then
     print("Years to analyze ("+syear+"-"+eyear+") fall outside temporal range of file ("+(syearF/100)+"-"+(eyearF/100)+"), exiting")
     exit
  end if
;-------------------------------------
; Many of the files have variables named the same as that indicated in the file name.
; Some however do not. Alter below.
;
  var = varT
  if (var.eq."zmta") then
     var = "ta"
  end if
  if (var.eq."monmaxpr") then
     var = "pr"
  end if
  if (var.eq."monmaxtasmax") then
     var = "tasmax"
  end if
  if (var.eq."monmintasmin") then
     var = "tasmin"
  end if  
;-------------------------------------
;  Read in filein, trim in time, remove annual cycle
;
  time0 := yyyymm_time(syearF/100,eyearF/100,"integer")  ; create time array matching time period specified in file name
  ind_s = ind(time0.eq.(syear*100+1))   ; identify analysis start index 
  ind_e = ind(time0.eq.eyear*100+12)   ; identify analysis end index 

  a = addfile(filein,"r")
  arr = a->$var$(ind_s:ind_e,:,:)  
  arr = rmMonAnnCycTLL(arr)
;-------------------------------------
;  Read in filein_tas, remove annual cycle
;
  a = addfile(filein_tas,"r")
  arr_tas = a->tas(ind_s:ind_e,:,:) 
  arr_tas = rmMonAnnCycTLL(arr_tas)
;-------------------------------------
;  Set up cosine of the latitude arrays for global average calculation, and compute global average of tas array (=G(t))          
;
  pi=4.*atan(1.0)
  rad=(pi/180.)
  coswgt    := cos(rad*arr_tas&lat)
  coswgt!0   = "lat"
  coswgt&lat = arr_tas&lat

  slat = -90.
  nlat = 90.
  global_aa := wgt_areaave(arr_tas(:,{slat:nlat},:),coswgt({slat:nlat}),1.0,0)     
;-------------------------------------
;  Write global average array to .nc file, call python script to calculate 10-yr butterworth low pass filter
;  and read smoothed timeseries back in.
;
  e = addfile("tmp.nc","c")
  e->aa = global_aa
  delete([/e,global_aa/])
  system("python butterworthfilter.aa.py")
  f = addfile("tmp2.nc","r")
  global_aaS = (/ f->aa /)
  system("rm tmp.nc tmp2.nc")
;-------------------------------------
;   Compute the linear regression of the smoothed global average tas timeseries against the input array anomalies
;
  global_rcS := regCoef_n(global_aaS,arr,0,0)
;-------------------------------------
;  Compute the linear trend of the smoothed G(t), and scale the regression coefficient by the linear trend of G(t)
;
  tttt := dtrend_msg_n(ispan(0,dimsizes(global_aaS)-1,1),global_aaS,True,True,0)  
  gtS := tttt@slope*dimsizes(global_aaS)   ; scale slope by number of timesteps
  global_rcS = (/ global_rcS*gtS /)     ; scale regression pattern by G(t) trend
  copy_VarMeta(arr(0,:,:),global_rcS)   ; copy metadata from arr array to global_rcS array)
;-------------------------------------
;  Plot resulting global_rcS array that contains the estimated forced response
;
  id = str_split(vT(2),".")      ; form output file name by parsing input file name
  if (temp(6).eq."Training") then
     rd = str_split(vT(5),".")
     if (temp(9).eq."CESM2") then   ; CESM2 filenames are slightly different
       runid = str_join(rd(:1),".")
     else
       runid = rd(0)
     end if
     ofile = id(0)+"_"+runid+"."+varT+".forced_response."+syear+"-"+eyear   ; output file name
  else
     ofile = id(0)+"."+varT+".forced_response."+syear+"-"+eyear   ; output file name
  end if
  wks = gsn_open_wks("ps",ofile)  ; open a workstation to hold the graphics

  res = True   ; set plot settings (resources)
  res@mpProjection = "WinkelTripel"
  res@mpGeophysicalLineColor = "gray42"
  res@mpGeophysicalLineThicknessF = 2.          
  res@mpPerimOn    = False
  res@mpGridLatSpacingF =  90            ; change latitude  line spacing
  res@mpGridLonSpacingF = 180.           ; change longitude line spacing
  res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
  res@mpGridAndLimbOn   = True           ; turn on lat/lon lines  
  res@mpFillOn = False
  res@mpCenterLonF = 210.
  res@mpOutlineOn = True  
  
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLineLabelsOn = False
  res@cnFillOn        = True
  res@cnLinesOn       = False
  res@lbLabelBarOn    = True
  res@pmLabelBarHeightF = 0.05
  res@pmLabelBarWidthF = 0.65
  res@lbLabelFontHeightF = 0.013

  res@gsnCenterStringOrthogonalPosF = 0.025
  res@gsnLeftStringParallelPosF = -.005
  res@gsnRightStringParallelPosF = 1.00
  res@gsnLeftStringOrthogonalPosF = -0.06
  res@gsnRightStringOrthogonalPosF = -0.05
  res@gsnLeftStringFontHeightF = 0.015
  res@gsnCenterStringFontHeightF = 0.019
  res@gsnRightStringFontHeightF = 0.015

  if (var.eq."psl") then    ; set contour levels depending on variable
     res@cnLevels := ispan(-500,500,50)
  end if
  if (var.eq."pr".or.var.eq."monmaxpr") then
     res@cnLevels := (/-6,-4,-2,-1,-0.5,-0.2,0,0.2,0.5,1,2,4,6/)
  end if
  if (var.eq."tas".or.var.eq."tos".or.var.eq."monmaxtasmax".or.var.eq."monmintasmin") then
     res@cnLevels := fspan(-2.5,2.5,21)
  end if
  if (var.eq."siconc") then
     res@cnLevels := ispan(-24,24,4)
  end if
  res@gsnCenterString = fn    ; set titles
  res@gsnLeftString = "G(t) Reg."
  res@gsnRightString = syear+"-"+eyear
  plot = gsn_csm_contour_map(wks,global_rcS,res)   ; plot global_rcS array
  delete(wks)
;-------------------------------------
; Use image magick to convert .ps file to png, then remove the .ps file.
;
  system("convert -density 196 -trim +repage -background white -border 10 -bordercolor white "+ofile+".ps "+ofile+".png")
  system("rm "+ofile+".ps")
  print(ofile+".png created!")
end

