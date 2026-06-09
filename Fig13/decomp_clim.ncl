; 计算气候态下气旋日降水和非气旋日降水

f = addfile("cyclone_flag.nc","r")
flag = f->flag
; printVarSummary(flag)

; AGCD降水
total_year = flag&year
nyear = dimsizes(total_year)
year_start = total_year(0)
year_end = total_year(nyear-1)
date_start = 817
date_end = 1215

nlat = 70
nlon = 89

cyclone_rain = new((/nlat,nlon/),"float")
cyclone_rain = 0.0

noncyclone_rain = new((/nlat,nlon/),"float")
noncyclone_rain = 0.0

D = 0. ;样本空间总日数

do year = year_start, year_end
    f = addfile("/g/data/if69/ls3248/data/AGCD/precip/daily/agcd_v1_precip_total_r05_daily_"+year+".nc","r")
    time = cd_calendar(f->time,2)
    ind_start = ind(time.eq.year*10000+date_start)
    ind_end = ind(time.eq.year*10000+date_end)
    precip = f->precip(ind_start:ind_end,:,:)
    delete(time)

    ;提取气旋日降水
    cyclone_ind = ind(flag({year},:).eq.1)
    cyclone_rain = cyclone_rain + dim_sum_n(precip(cyclone_ind,:,:), 0)

    ;提取非气旋日降水
    noncyclone_ind = ind(flag({year},:).eq.0)
    noncyclone_rain = noncyclone_rain + dim_sum_n(precip(noncyclone_ind,:,:), 0)

    D = D + dimsizes(cyclone_ind) + dimsizes(noncyclone_ind)

    delete([/cyclone_ind,noncyclone_ind/])
    print(year+" done")
end do

pc_clim = cyclone_rain/D
pnc_clim = noncyclone_rain/D

pc_clim!0 = "lat"
pc_clim!1 = "lon"
pc_clim&lat = f->lat
pc_clim&lon = f->lon  
copy_VarCoords(pc_clim, pnc_clim)

;-----output the Comp_avg_hgt
out_file = "decomp_clim.nc"

if (fileexists(out_file)) then
    system("rm -f "+out_file)
end if
f_out = addfile(out_file,"c")
f_out->pc_clim = pc_clim
f_out->pnc_clim = pnc_clim

;plot test
term = new((/2,nlat,nlon/),"float")
term(0,:,:) = pc_clim
term(1,:,:) = pnc_clim

latMin = -45
latMax = -10
lonMin = 105
lonMax = 155

wks = gsn_open_wks("pdf","decomp_clim")  ; open a workstation
	res = True
	res@gsnDraw = False
	res@gsnFrame = False
	res@gsnAddCyclic = False
	res@gsnMaximize = False
	res@cnFillOn = True
	res@cnLinesOn = False
	res@cnLineLabelsOn = False
	res@cnInfoLabelOn = False
    res@lbLabelBarOn         = False
	res@tmXTOn = False
	res@tmYROn = False

	res@mpMinLatF = latMin 
    res@mpMaxLatF = latMax
    res@mpMinLonF = lonMin   
    res@mpMaxLonF = lonMax
    res@mpCenterLonF = (lonMin+lonMax)/2.
	res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = fspan(0, 2, 11)
    res@cnFillPalette = "cmocean_deep"
	res@tmYLLabelFontHeightF = 0.03
	res@tmXBLabelFontHeightF = 0.03
    res@mpFillOn = True
    res@mpDataSetName = "Earth..4"
    res@mpOutlineBoundarySets = "AllBoundaries"
	res@mpProvincialLineColor  = "grey"
	res@mpProvincialLineThicknessF = 1.5
    res@cnFillDrawOrder = "PreDraw"
    res@mpOceanFillColor = "white" 
    res@mpFillDrawOrder = "Draw"
    res@mpLandFillColor = "transparent"
    res@mpDataBaseVersion  = "MediumRes"
    res@mpOutlineDrawOrder = "Draw"
	res@mpShapeMode = "FreeAspect"
	res@gsnMajorLonSpacing = 10

	centerstring = (/"cyclone day","non-cyclone day"/)
	plot = new(2, graphic)
    do i = 0, 1
        res@gsnCenterString = centerstring(i)
		plot(i) = gsn_csm_contour_map_ce(wks,term(i,:,:),res)
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	; resP@pmLabelBarWidthF = 0.035
	; resP@pmLabelBarHeightF = 0.6
	; resP@lbOrientation = "Vertical"
	resP@pmLabelBarOrthogonalPosF = 0.01
	resP@lbLabelFontHeightF = 0.012
	gsn_panel(wks,plot,(/1,2/),resP)  
