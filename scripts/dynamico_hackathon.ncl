/; This code is part of the dynamico set of dynamical adjustment codes
---
--- Program dynamico_hackathon.ncl
--- This code performs dynamical adjustment as described in Deser et al. 2016 J.Climate.
--- The code is an example for the August 2023 ForceSMIP Hackathon at ETH Zurich.
---
--- The code has to be used with the ForceSMIP training dataset
--- 
--- The code reads variables var (e.g tas or pr) and varp (e.g psl) from one member of the model training dataset, 
--- detrend variables var and/or varp using a parametric or non parametric
--- method and then use the constructed analogue method to carry out dynamical adjustment. The code outputs the var dynamic component,
--- the reconstructed varp as well as all iterations of the var dynamic component (can be used to assess uncertainties).
--- In Evaluation mode, it also outputs the ForceSMIP "forced component" which is simply the linear trend of the residual annual means 
--- (the latter being obtained by subtracting the dynamic component from the original var array). 
---
--- Command-line parameters:
--- General: MODE, MODEL, EXPnm, var, varp, member, tStart, tEnd, ANOM, ANOMP, FMTnm
--- *** MODE: run in Training or Evaluation mode
--- *** MODEL: one of the ForceSMIP models (e.g MIROC-ES2L)
--- *** EXPnm: the name of the simulation for the above model (e.g historical_ssp245)
--- *** tStart, tEnd: beginning and ending of the ForceSMIP period (used to define ForceSMIP filenames) 
--- *** member: the full member "name" (e.g r1i1p1f2 or r1031.002i1p1f1 for CESM2)
--- *** var: the variable to be dynamically adjusted (e.g tas or pr)
--- *** varp: the atmospheric circulation variable (e.g psl)
--- *** ANOM: use anomalies instead of raw values for var (climatology is estimated over the full period)
--- *** ANOMP: use anomalies instead of raw values for varp (climatology is estimated over the full period)
--- *** region: geographic domain to apply dynamical adjustment. Currently: GLO (global), NH (northern hemisphere) NA (North America) and EU (Europe)
---             Other domains can easily be added by users (see section 1. of the code below)
---
--- Dynamical adjustment method: N_a, N_b, niter
--- *** N_a: maximum number of a specific month (if N years present, then set N_a to N-1. For ForceSMIP Training, N_a = 142)
--- *** N_b: number of months to be sub-sampled for the reconstruction (note that N_b must be < N_a). 70-80% of N_a is a good starting value.
--- *** niter: number of iterations for the dynamical adjustment algorithm (usually between 30 and 100). Use a small number first -3- to test the code.
--- 
--- Detrending: DTRD, DTRDP, nl, nsmo, knots
--- *** DTRD: detrending options for var (none, linear, quadratic, cubic, spline, loess, gsat_loess). Note that the use of loess and gsat_loess requires an extra package (NCSTAT)
--- *** DTRDP: detrending options for varp (none, linear, quadratic, cubic, spline)
--- *** nl, smo: parameters of the loess routine
--- *** knots: number of knots for the cubic spline (for training period 1880-2022, knots value should be between 6 and 8)
---
--- ASSUMPTIONS: 1. All data files (psl, tas, pr)  must have the same order for lon (0->360)
---                 as the code assumes it is the case (use cdo or nco tricks to do it beforehand if needed)
---                 For latitude, the code assumes that the latitudes go as (-90->90). If not the case, use cdo invertlat beforehand.
---              2. The code assumes that all data files have correct time coordinates (can be read by cd_calendar)
---              3. The code assumes that the data have no missing values so no infilling is performed 
---
--- Calling line examples: 
--- Training mode: ncl -Q var=\"tas\" varp=\"psl\" MODE=\"Training\" MODEL=\"CanESM5\" EXPnm=\"historical_ssp585\" tStart=188001 tEnd=202212 member=\"r1i1p2f1\"  N_a=142 N_b=100 niter=30  region=\"EU\" DTRD=\"quadratic\" DTRDP=\"none\" nl=45 nsmo=5 knots=7 ANOM=False  ANOMP=False FMTnm=\"NetCDF4Classic\" dynamico_hackathon.ncl
---
--- Evaluation mode: ncl -Q var=\"tas\" varp=\"psl\" MODE=\"Evaluation\" TIER=\"Tier1\" MODEL=\"1H\" EXPnm=\"Unknown\" tStart=195001 tEnd=202212 member=\"Unknown\"  N_a=72 N_b=60 niter=30  region=\"GLO\" DTRD=\"quadratic\" DTRDP=\"none\" nl=45 nsmo=5 knots=7 ANOM=False  ANOMP=False FMTnm=\"NetCDF4Classic\" dynamico_hackathon.ncl
---
--- To start, the user has to put this code in the ForceSMIP root directory (with the Training, Evaluation-Tier1,  Evaluation-Tier2, Evaluation-Tier3 sub-directories).
--- Then the user has to create the work and results directories (see Input section just below). You are then ready to go !
---
--- Author: Laurent Terray 08/20/2023
----------------------------------------------------------------------------------------------------------------------
;/
begin
;---------------------------------------------------------------------------------------------------------------------
;--- Input: ROOT is the directory with the ForceSMIP datasets
    ROOT = "/scratch/globc/terray/ForceSMIP/" ; root directory
    DIRO = ROOT+"results/"  ; directory to save output file
    WORK = ROOT+"work/"     ; temporary work directory
;--- End of Input ----------------------------------------------------------------------------------------------------
;
;--- Section 1: here the user can add new geographical domains
;
;--- Regions: if the longitude domain crosses 0 then flip must equal to "yes"  
    if (region .eq. "EU") then
      flip = "yes"
      minlat = 30
      maxlat = 80
      minlon = -20
      maxlon = 50
    end if
    if (region .eq. "GLO") then
      flip = "no"
      minlat = -90
      maxlat = 90
      minlon = 0
      maxlon = 360
    end if
    if (region .eq. "NA") then
      flip = "no"
      minlat = 20
      maxlat = 80
      minlon = 190
      maxlon = 320
    end if
    if (region .eq. "NH") then
      flip = "no"
      minlat = 0
      maxlat = 90
      minlon = 0
      maxlon = 360
    end if
;---------------------------------------------------------------------------------------------------------------------
;---
;--- Section 2: read data and perform basic checks
    period = tStart+"-"+tEnd
;--- Load data from historical simulation
    if (MODE .eq. "Training") then
      rprefix=ROOT+MODE+"/Amon/"+varp+"/"+MODEL+"/"
      f_rpsl  = addfile(rprefix+varp+"_mon_"+MODEL+"_"+EXPnm+"_"+member+"."+period+".nc","r")
      rprefix=ROOT+MODE+"/Amon/"+var+"/"+MODEL+"/"
      f_rtas  = addfile(rprefix+var+"_mon_"+MODEL+"_"+EXPnm+"_"+member+"."+period+".nc","r")
    end if
    if (MODE .eq. "Evaluation") then
      rprefix=ROOT+MODE+"-"+TIER+"/Amon/"+varp+"/"
      f_rpsl  = addfile(rprefix+varp+"_mon_"+MODEL+"."+period+".nc","r")
      rprefix=ROOT+MODE+"-"+TIER+"/Amon/"+var+"/"
      f_rtas  = addfile(rprefix+var+"_mon_"+MODEL+"."+period+".nc","r")
    end if
;---
;--- Loading variables on training period
;---
    print("   Loading circulation variable "+varp+" and the other variable "+var+" on the full period")
;--- Read circulation variable and adjusted variable (use float instead of double)
    tmp  = dble2flt(f_rpsl->$varp$)
;--- Define scale factors for units
    if (varp .eq. "psl") then
      scale = .01               ; scale psl (back to hPa)
    else
      scale = 1.
    end if
    tmp = (/ tmp * scale/)
    tmp_t = dble2flt(f_rtas->$var$)
    unit_nm = tmp_t@units
;--- time stuff
    time = f_rpsl->time
    cldr = cd_calendar(time, -1)
    ntime = dimsizes(cldr)
    nyrs = ntime/12
    print("   File initial date = "+cldr(0)+" - Final date = "+cldr(ntime-1)+" - Total number of years = "+nyrs)
    delete([/ cldr/])
;--- Check if longitudes and latitudes are correctly ranked
    lon_g = f_rpsl->lon
    lat_g = f_rpsl->lat
    nlat = dimsizes(lat_g)
    if (min(lon_g) .lt. 0) then
      print("   Some longitudes are negative : "+min(lon_g) )
      print("   Longitudes must be within 0 - 360; The code stops")
      exit
    end if
    if (lat_g(0) .gt. lat_g(nlat-1)) then
      print("   The latitude array goes from North to South")
      print("   The code requires the opposite; the code stops")
      exit
    end if
    delete([/lat_g, lon_g/])
;--- Checks
    print("   Min/Max of circulation variable "+varp+" = "+min(tmp)+" - "+max(tmp) )
    print("   Min/Max of the other variable "+var+" = "+min(tmp_t)+" - "+max(tmp_t) )
    if (flip .eq. "yes") then
      tmp = lonFlip(tmp)
      tmp_t = lonFlip(tmp_t)
    end if
;--- Get anomalies if needed
    if (ANOMP) then
      tmp_clm     = clmMonTLL(tmp)
      tmp  = (/ calcMonAnomTLL(tmp,tmp_clm) /)
      delete(tmp_clm)
      print("   Min/Max "+varp+" anoomalies:")
      printMinMax(tmp, True)
    end if
    if (ANOM) then
      tmp_clm     = clmMonTLL(tmp_t)
      tmp_t  = (/ calcMonAnomTLL(tmp_t,tmp_clm) /)
      delete(tmp_clm)
      print("   Min/Max "+var+" anoomalies:")
      printMinMax(tmp_t, True)
    end if
;--- If Evaluation, keep a copy of the original data to get the residual at the end
    if (MODE .eq. "Evaluation") then
      tmp_t_sav = tmp_t
    end if
;---------------------------------------------------------------------------------------------------------------------
;---
;--- Section 3: detrending step
;--- 
    print(DTRD+" detrending of "+var+" performed for each month")
    if (DTRD .eq. "linear") then
      do imo = 0,11
        tmp_t(imo::12,:,:) = (/ dtrend_msg_n(ispan(1,nyrs,1),tmp_t(imo::12,:,:), False, False, 0)  /)
      end do
    elseif (DTRD .eq. "quadratic") then
      do imo = 0,11
        tmp_t(imo::12,:,:) = (/ dtrend_quadratic_msg_n(tmp_t(imo::12,:,:), False, False, 0)  /)
      end do
    elseif (DTRD .eq. "cubic") then
      dt = ispan(1,nyrs,1)
      ndeg = 4  ;--- number of coefficients searched
      wgts = 1. ;--- equal weights
      do imo = 0,11
        mean = dim_avg_n_Wrap(tmp_t(imo::12,:,:),0)
        coef = lspoly_n(dt, tmp_t(lat|:,lon|:,time|imo::12), wgts, ndeg, 2)
        poly3 = tmp_t(imo::12,:,:)
        do it = 0, nyrs -1
          poly3(it,:,:) = (/coef(:,:,0) + coef(:,:,1) * dt(it) + coef(:,:,2) * dt(it)^2 + coef(:,:,3) * dt(it)^3 /)
        end do 
        cubic = tmp_t(imo::12,:,:)
        cubic = (/ tmp_t(imo::12,:,:) - poly3 + conform_dims(dimsizes(tmp_t(imo::12,:,:)), mean, (/ 1,2 /) )/)
        tmp_t(imo::12,:,:) = (/ cubic /)
        delete ([/ cubic, poly3, coef, mean /])
      end do
    elseif (DTRD .eq. "spline") then
      print("   spline detrending of "+var+" is performed")
      xi = fspan(1, nyrs,nyrs)
      yo = tmp_t(lat|:,lon|:,time|0::12)
      do imo = 0,11
        yz = tmp_t(lat|:,lon|:,time|imo::12)
        yo = (/ csa1s(xi,yz,knots,xi) /)
        mean = dim_avg_n_Wrap(tmp_t(imo::12,:,:), 0)
        tmp_t(imo::12,:,:) = (/ tmp_t(imo::12,:,:) - yo(time|:,lat|:,lon|:) /)
        tmp_t(imo::12,:,:) = (/ tmp_t(imo::12,:,:) + conform_dims(dimsizes(tmp_t(imo::12,:,:)), mean,(/1,2/) ) /)
        delete([/ yz, mean /])
      end do
      delete([/ yo, xi /])
    elseif (DTRD .eq. "loess") then
;--- Get specific string for file naming
      datestr = systemfunc("date")
      datestr = str_sub_str(datestr, " ", "_")
      datestr = str_sub_str(datestr, ":", "_")
;--- Detrend each month separately (to account for seasonality of trends
      finp_nm = WORK+var+"_inp_"+member+datestr+".nc"
      fout_nm = WORK+var+"_out_"+member+datestr+".nc"
      v_nm = var+"_"+"residual"
      w_nm = var+"_"+"trend"
      system("rm -rf " + finp_nm)
      system("rm -rf " + fout_nm)
      do imo = 0, 11
        tmp_mo = tmp_t(imo::12,:,:)
        print(" Min/Max "+var+" before loess filter for month:"+(imo+1)+" = "+min(tmp_mo)+" - "+max(tmp_mo) ) 
        dim_t = dimsizes(tmp_mo)
        ntim = dim_t(0)
        finp = addfile(finp_nm,"c")
        filedimdef(finp,"time",-1,True)
        finp->$var$ = tmp_mo
        cmd = "comp_trend_3d -f="+finp_nm+" -v="+var+" -nt="+nl+" -smt="+nsmo+" -a=residual -robust -o="+fout_nm
        system(cmd)
        fout = addfile(fout_nm,"r")
        trdd = fout->$w_nm$
        resi = fout->$v_nm$
        resi@_FillValue = tmp_t@_FillValue
        print(" 1. Min/Max residual "+var+" after loess filter for month:"+(imo+1)+" = "+min(resi)+" - "+max(resi) )
;--- Add trend value at initial date to the residual
        tmp_gm = trdd(0,:,:)
        resi = (/ resi + conform_dims(dimsizes(resi), tmp_gm, (/1, 2/) ) /)
        print(" 2. Min/Max full "+var+" after loess filter for month:"+(imo+1)+" = "+min(resi)+" - "+max(resi) )
        tmp_t(imo::12,:,:) = (/ resi /)
        system("rm -rf " + finp_nm)
        system("rm -rf " + fout_nm)
        delete([/tmp_mo, dim_t, ntim, finp, fout, resi, cmd, trdd, tmp_gm/])
      end do
    elseif (DTRD .eq. "gsat_loess") then
      ntl = nl * 12                        ;--- Number of months for the the loess filter length
      smt = nsmo * 12                     ;--- Number of months for the additional smoothing
;--- Get global mean temperature   
      lat_g = dble2flt(f_rtas->lat)
      wgt = latRegWgt(lat_g, "float", 0)
      if (var .eq. "tas") then
        gsat = wgt_areaave_Wrap(tmp_t, wgt, 1., 0)
      else
        print(" to use gsat_loess for "+var+" one first need to read the tas variable to calculate gsat")
        tprefix=ROOT+MODE+"/Amon/tas/"+MODEL+"/"
        ftas  = addfile(tprefix+"tas_mon_"+MODEL+"_"+EXPnm+"_"+member+"."+period+".nc","r")
        tas = dble2flt(ftas->tas)
        gsat = wgt_areaave_Wrap(tas, wgt, 1., 0)
        delete(tas)
      end if
;--- Smooth GSAT Loess estimate
      resi_nm = var + "_residual"
      trd_nm = var + "_trend"
      finp_nm = WORK+"Input_GSAT_1m_"+MODEL+"_"+member+".nc"
      fout_nm = WORK+"Output_GSAT_1m_"+MODEL+"_"+member+".nc"
      system("rm -rf " + finp_nm)
      system("rm -rf " + fout_nm)
      finp = addfile(finp_nm,"c")
      filedimdef(finp,"time",-1,True)
      finp->$var$ = gsat
      cmd = "comp_trend_1d -f="+finp_nm+" -v="+var+" -nt="+ntl+" -smt="+smt+" -a=residual -robust -o="+fout_nm
      system(cmd)
      fout = addfile(fout_nm,"r")
      gsat_flt = fout->$trd_nm$
;--- Get delta GSAT 
      nmonth = dimsizes(gsat_flt)
      gwll = gsat_flt(nmonth-1) - gsat_flt(0)
      gwll = decimalPlaces(gwll, 2, True)
      print(" Delta in GSAT for the period "+period+" = "+gwll)
;--- Perform the regression on GSAT
      index_nm = "gsat"
      fvar_nm = WORK+"fichvar_"+MODEL+"_"+member+".nc"
      find_nm = WORK+"fichind_"+MODEL+"_"+member+".nc"
      freg_nm = WORK+"fichreg_"+MODEL+"_"+member+".nc"
      resid_nm = var+"_"+index_nm+"_resid"
      system("rm -rf " + fvar_nm)
      system("rm -rf " + find_nm)
      system("rm -rf " + freg_nm)
      fvar = addfile(fvar_nm,"c")
      filedimdef(fvar,"time",-1,True)
      find = addfile(find_nm,"c")
      filedimdef(find,"time",-1,True)
      find->$index_nm$ = gsat_flt
      fvar->$var$ = tmp_t
      cmd = "comp_reg_3d -f="+fvar_nm+" -v="+var+" -p=1 -fi="+find_nm+" -vi="+index_nm+" -pi=1,1 -a=residual -add_mean -o="+freg_nm
      system(cmd)
      freg = addfile(freg_nm,"r")
      resi = freg->$resid_nm$
      printVarSummary(resi)
      print(" Min/max "+var+" after regression = "+min(resi)+" - "+max(resi) )
      tmp_t = (/ resi /)
      system("rm -rf " + finp_nm)
      system("rm -rf " + fout_nm)
      system("rm -rf " + fvar_nm)
      system("rm -rf " + find_nm)
      system("rm -rf " + freg_nm)
      delete([/ resi, fvar, find, freg/])
    elseif (DTRD .eq. "none") then
      print("   No detrending of "+var)
    else
      print("   detrending not implemented; Stop the code")
      exit
    end if
    if (DTRD .ne. "none") then
      print("   Min/Max after detrending of variable "+var+" = "+min(tmp_t)+" - "+max(tmp_t) )
    end if
;--- Detrending of circulation variable
    print(DTRDP+"  detrending of "+varp+" performed for each month")
    if (DTRDP .eq. "linear") then
      nyrs = ntime/12
      do imo = 0,11
        tmp(imo::12,:,:) = (/ dtrend_msg_n(ispan(1,nyrs,1),tmp(imo::12,:,:), False, False, 0)  /)
      end do
    elseif (DTRDP .eq. "quadratic") then
      do imo = 0,11
        tmp(imo::12,:,:) = (/ dtrend_quadratic_msg_n(tmp(imo::12,:,:), False, False, 0)  /)
      end do
    elseif (DTRDP .eq. "cubic") then
      dt = ispan(1,nyrs,1)
      ndeg = 4  ;--- number of coefficients searched
      wgts = 1. ;--- equal weights
      do imo = 0,11
        coef = lspoly_n(dt, tmp(lat|:,lon|:,time|imo::12), wgts, ndeg, 2)
        poly3 = tmp(imo::12,:,:)
        do it = 0, nyrs -1
          poly3(it,:,:) = (/coef(:,:,0) + coef(:,:,1) * dt(it) + coef(:,:,2) * dt(it)^2 + coef(:,:,3) * dt(it)^3 /)
        end do 
        cubic = tmp(imo::12,:,:)
        cubic = (/ tmp(imo::12,:,:) - poly3 /)
        tmp(imo::12,:,:) = (/ cubic /)
        delete ([/ cubic, poly3, coef /])
      end do
    elseif (DTRDP .eq. "spline") then
      xi = fspan(1, nyrs,nyrs)
      yo = tmp(lat|:,lon|:,time|0::12)
      do imo = 0,11
        yz = tmp(lat|:,lon|:,time|imo::12)
        yo = (/ csa1s(xi,yz,knots,xi) /)
        tmp_clm = dim_avg_n_Wrap(tmp(imo::12,:,:), 0)
        tmp(imo::12,:,:) = (/ tmp(imo::12,:,:) - yo(time|:,lat|:,lon|:) /)
        tmp(imo::12,:,:) = (/ tmp(imo::12,:,:) + conform_dims(dimsizes(tmp(imo::12,:,:)), tmp_clm,(/1,2/) ) /)
        delete([/ yz, tmp_clm /])
      end do
      delete([/ yo, xi /])
    elseif (DTRDP .eq. "none")
      print(" ")
    else
      print("   detrending not implemented; Stop the code")
      exit
    end if
    if (DTRDP .ne. "none") then
      print("   Min/Max after detrending of circulation variable "+varp+" = "+min(tmp)+" - "+max(tmp) )
    end if
;--- Get final arrays before dynamical adjustment
    rpsl_all = tmp(time|:, {lat|minlat:maxlat},{lon|minlon:maxlon})
    rtas_all = tmp_t(time|:, {lat|minlat:maxlat},{lon|minlon:maxlon})
    rpsl = tmp(time|12:, {lat|minlat:maxlat},{lon|minlon:maxlon})
    rtas = tmp_t(time|12:, {lat|minlat:maxlat},{lon|minlon:maxlon})
    delete([/ tmp, tmp_t, f_rtas, f_rpsl/])
;---------------------------------------------------------------------------------------------------------------------
;---
;--- Section 4: Pre-processing before dynamical adjustment
;---
;--- Getting latitudes & longitudes
    lat              = rpsl&lat
    lon              = rpsl&lon
    nlat             = dimsizes(lat)
    nlon             = dimsizes(lon)
;--- Create arrays to hold dynamically adjusted variable and reconstructed circulation analogues
    tas_adj=conform_dims( (/niter,ntime,nlat,nlon/), rtas_all, (/1,2,3/) )
    itera=ispan(1,niter,1)
    itera!0 = "iteration"
    tas_adj!0 = "iteration"
    tas_adj&iteration = itera
    tas_adj!1 = "time"
    tas_adj&time = rtas_all&time
    tas_adj!2 = "lat"
    tas_adj&lat = lat
    tas_adj!3 = "lon"
    tas_adj&lon = lon
;--- Circulation variable now
    psl_rec = conform_dims( (/niter,ntime,nlat,nlon/), rpsl_all, (/1,2,3/) )
    psl_rec!0 = "iteration"
    psl_rec&iteration = itera
    psl_rec!1 = "time"
    psl_rec&time = rpsl_all&time
    psl_rec!2 = "lat"
    psl_rec&lat = lat
    psl_rec!3 = "lon"
    psl_rec&lon = lon
;--- Add attributes
    today = systemfunc("date")
    tas_adj@long_name = rtas_all@long_name+" Dynamical contribution for all iterations"
    psl_rec@long_name = "Constructed "+rpsl_all@long_name+" analog for all iterations"
    psl_rec@description = " Dynamico run performed on "+today(0)+" with algorithm parameters: N_a = "+N_a+" , N_b = "+N_b+" , niter = "+niter+" ; Model: "+MODEL+" Member: "+member+" - Historical simulation: "+EXPnm+" years: "+tStart+"-"+tEnd
    tas_adj@description = psl_rec@description
    ntt = nyrs -1
;---------------------------------------------------------------------------------------------------------------------
;---
;--- Section 5: The leave one out method 
;---    fpsl is the target array, rpsl the training analogue array
;---
    print("  Constructed analogue algorithm: time loop begins")
    do ii = 0, ntime - 1
      iyear = ii/12
      iyearn = tStart/100 + iyear
      isea = ii - (iyear*12)
      print("   time step = "+ii+" - Year: "+iyearn+" - Month: "+(isea+1) )
;--- get the array to be reconstructed
      fpsl = rpsl_all(time|ii, lat|:, lon|:)
      etas = rtas_all(time|ii, lat|:, lon|: )
;--- get the training sample and deal with boundaries
      if (iyear .eq. 0) then      ; target is in first year
        rpsl = (/ rpsl_all(time|12:, lat|:, lon|:) /)
        rtas = (/ rtas_all(time|12:, lat|:, lon|:) /)
      elseif (iyear .eq. ntt) then    ; target is in last year
        nend = ntime - 13
        rpsl = (/ rpsl_all(time|:nend, lat|:, lon|:) /)
        rtas = (/ rtas_all(time|:nend, lat|:, lon|:) /)
      else  ; normal month
        imb = iyear*12 - 1
        ima = iyear*12 + 12
        rpsl(:imb,:,:) =  (/ rpsl_all(time|:imb, lat|:, lon|:) /)
        rpsl(imb+1:,:,:) =  (/ rpsl_all(time|ima:, lat|:, lon|:) /)
        rtas(:imb,:,:) =  (/ rtas_all(time|:imb, lat|:, lon|:) /)
        rtas(imb+1:,:,:) =  (/ rtas_all(time|ima:, lat|:, lon|:) /)
      end if
;--- Select only the current month analogues
      dpsl = rpsl(isea::12,:,:)
      dtas = rtas(isea::12,:,:)
/;-----------------------------------------------------------
 Constructed analog step using Moore-Penrose pseudo-inverse
 and perform loop on iteration of analog steps
;/
      fpsl_m = onedtond(ndtooned(fpsl),(/1,nlat*nlon/))  ;--- Put target field in vector shape
      rec_psl_2d = new((/niter,nlat,nlon/), float, rpsl@_FillValue) 
      rec_tas_2d = new((/niter,nlat,nlon/), float, rtas@_FillValue)
      do ia = 0,niter-1         ; Start of iteration loop
;--- Select randomly N_b analogs among the N_a ones
        nba = generate_unique_indices(N_a)
        iadr = nba(0:N_b-1)
        A_mat = onedtond(ndtooned(dpsl(iadr,:,:)),(/N_b,nlat*nlon/))
        dimss = dimsizes(A_mat)
        nca = dimss(0)
        nra = dimss(1)
        s = new((/nca,nca/),float)
        s = (/ 0. /)
        u   = new ( (/nca,nra/) , typeof(A_mat) )  ;--- pre-allocate space 
        v   = new ( (/nca,nca/) , typeof(A_mat) )
        sgv  = svd_lapack (A_mat, "S" , "S", 1,  u, v)  ;--- SVD step - if 0 v is here v_T
        delete([/ A_mat, dimss /])
        do im = 0, nca - 1
          s(im,im) = (/ 1./sgv(im) /)
        end do
;--- Perform pseudo-inverse to get linear coefficients and apply to analogues
        x = fpsl_m(0,:) # transpose(u) # s # v       ;--- x has dimension N_b
        delete([/ s, u, v, sgv, nba/])
        rec_tas = x # onedtond(ndtooned(dtas(iadr,:,:)),(/N_b,nlat*nlon/))
        rec_psl = x # onedtond(ndtooned(dpsl(iadr,:,:)),(/N_b,nlat*nlon/))
        delete([/ x, iadr/])
        rec_tas_2d(ia,:,:) = onedtond(rec_tas,(/nlat,nlon/))
        rec_psl_2d(ia,:,:) = onedtond(rec_psl,(/nlat,nlon/))
        delete([/ rec_tas, rec_psl/])
      end do               ;--- end of loop on number of iterations
      tas_adj(:,ii,:,:) = (/rec_tas_2d /)
      psl_rec(:,ii,:,:) = (/rec_psl_2d /)
      delete([/ rec_tas_2d, rec_psl_2d, fpsl_m, fpsl, dtas, dpsl/])
    end do                 ; end of time loop
;--- 
    print("   CPU time: " + get_cpu_time())
    print("   CPU time: " + get_cpu_time())
    print("   Analog reconstruction done for model "+MODEL)
;
;--- Average over niter iterations
    tas_adj_avg = dim_avg_n_Wrap(tas_adj,0)
    psl_rec_avg = dim_avg_n_Wrap(psl_rec,0)
    print("   Write Output files in folder "+DIRO)
;--- If Evaluation mode, get the residual and calculate the linear trend
    if (MODE .eq. "Evaluation") then
      forced_trend = tmp_t_sav(0,:,:)
      residual = tmp_t_sav
      residual = (/ tmp_t_sav - tas_adj_avg /)
;--- get annual means
      nbm = 12
      inm = 5
      f_rm = runave_n_Wrap(residual,nbm,0,0)
      residual := f_rm(inm::12,:,:)
      dims = dimsizes(residual)
      years = ispan(1, nyrs, 1)
      printVarSummary(residual)
;--- get linear trend
      trd = dtrend_msg_n (years,residual,False,True,0)
      forced_trend = (/ onedtond(trd@slope, (/dims(1),dims(2)/) ) * tofloat(nyrs) /)
      forced_trend@long_name = "Forced component of "+var+" for model "+MODEL
      forced_trend@standard_name = "Forced component"
      forced_trend@units = unit_nm+" per "+nyrs+" years"
      delete(forced_trend@time)
      printVarSummary(forced_trend)
    end if
;---------------------------------------------------------------------------------------------------------------------
/;
-------------------------
---- Section 6.  Write results
-------------------------
;/
;--- Write output data
    var_a = var + "_avg"
    var_i = var + "_all"
    varp_a = varp + "_avg"
    varp_i = varp + "_all"
    var_f = "forced_component"
;--- Opening output files
    setfileoption("nc", "Format",  FMTnm)
    if (MODE .eq. "Training") then
      filo = var+"_"+varp+"_"+MODEL+"_member_"+member+"_simulation_"+EXPnm+"_period_"+period+"_parameters_"+N_a+"_"+N_b+"_"+niter+"_region_"+region+"_detrend_"+DTRD+"-"+DTRDP+".nc"
    elseif (MODE .eq. "Evaluation") then
      filo = var+"_"+varp+"_model_"+MODEL+"_period_"+period+"_parameters_"+N_a+"_"+N_b+"_"+niter+"_region_"+region+"_detrend_"+DTRD+"-"+DTRDP+".nc"
    end if
    system("/bin/rm -f " + DIRO + filo)
    fout  = addfile (DIRO + filo, "c")
    setfileoption(fout,"DefineMode",True)
;--- file for adjusted variable first
    dimNames = (/"iteration","lat", "lon", "time"/)  
    dimSizes = (/ niter, nlat,  nlon, ntime/) 
    dimUnlim = (/ False, False, False, True/)   
    filedimdef(fout,dimNames,dimSizes,dimUnlim)
    filevardef(fout, "iteration",typeof(itera),getvardims(itera))
    filevardef(fout, "lat"  ,typeof(lat),getvardims(lat))
    filevardef(fout, "lon"  ,typeof(lon),getvardims(lon))
    filevardef(fout, "time"  ,typeof(time),getvardims(time))
    filevardef(fout, var_i    ,typeof(tas_adj)  ,getvardims(tas_adj))
    filevardef(fout, var_a    ,typeof(tas_adj_avg)  ,getvardims(tas_adj_avg))
    filevardef(fout, varp_i    ,typeof(psl_rec)  ,getvardims(psl_rec))
    filevardef(fout, varp_a    ,typeof(psl_rec_avg)  ,getvardims(psl_rec_avg))
    if (MODE .eq. "Evaluation") then
      filevardef(fout, var_f    ,typeof(forced_trend)  ,getvardims(forced_trend))
      filevarattdef(fout, var_f  ,forced_trend)
    end if
    filevarattdef(fout,"lat"  ,lat)
    filevarattdef(fout,"lon"  ,lon)
    filevarattdef(fout,"time"  ,time)
    filevarattdef(fout,var_a  ,tas_adj_avg)
    filevarattdef(fout,var_i  ,tas_adj)
    filevarattdef(fout,varp_a  ,psl_rec_avg)
    filevarattdef(fout,varp_i  ,psl_rec)
    setfileoption(fout,"DefineMode",False)
    fout->iteration = (/itera/)
    fout->lat    = (/lat/)
    fout->lon    = (/lon/)
    fout->time   = (/time/)
    fout->$var_a$     = (/tas_adj_avg/)
    fout->$var_i$     = (/tas_adj/)
    if (MODE .eq. "Evaluation") then
      fout->$var_f$     = (/forced_trend/)
    end if
    fout->$varp_a$     = (/psl_rec_avg/)
    fout->$varp_i$     = (/psl_rec/)
    delete(fout)
;********************************************************************************************
end
