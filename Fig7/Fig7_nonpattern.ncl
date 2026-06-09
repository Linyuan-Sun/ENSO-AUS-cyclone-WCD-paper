
sphere = "SH"
latMax = 0
latMin = -90;-60
lonMax = 360;180
lonMin = 0;90
siglvl = 0.5 ;sig :0/1

season = "SON"
pltTpye = "pdf"

pattern = "nonpattern" ;"pattern","nonpattern"

if pattern.eq."pattern" then
	pltName = "Fig7_"+pattern
	f = addfile("cyclone_pattern_"+season+".nc", "r") ;total days in each month
	f_EP = addfile("./bootstrap/bootstrap_sig_EPElNino_pattern_SON.nc", "r") 
	f_CP = addfile("./bootstrap/bootstrap_sig_CPElNino_pattern_SON.nc", "r") 
	f_La = addfile("./bootstrap/bootstrap_sig_LaNina_pattern_SON.nc", "r") 
elseif pattern.eq."nonpattern"
	pltName = "Fig7_"+pattern
	f = addfile("cyclone_nonpattern_"+season+".nc", "r") ;total days in each month
	f_EP = addfile("./bootstrap/bootstrap_sig_EPElNino_nonpattern_SON.nc", "r") 
	f_CP = addfile("./bootstrap/bootstrap_sig_CPElNino_nonpattern_SON.nc", "r") 
	f_La = addfile("./bootstrap/bootstrap_sig_LaNina_nonpattern_SON.nc", "r") 
end if

EPElNino = f->freq_EPElNino({latMin:latMax},:)
CPElNino = f->freq_CPElNino({latMin:latMax},:)
LaNina = f->freq_LaNina({latMin:latMax},:)
freq_clim = f->freq_clim({latMin:latMax},:)
EPElNino_anom = EPElNino - freq_clim
CPElNino_anom = CPElNino - freq_clim
LaNina_anom = LaNina - freq_clim
copy_VarCoords(EPElNino,EPElNino_anom)
copy_VarCoords(CPElNino,CPElNino_anom)
copy_VarCoords(LaNina,LaNina_anom)

EPElNino_prob = f_EP->sig_mask({latMin:latMax},:)
CPElNino_prob = f_CP->sig_mask({latMin:latMax},:)
LaNina_prob = f_La->sig_mask({latMin:latMax},:)

    ;-----comp
	nlat = dimsizes(EPElNino_anom(:,0))
	nlon = dimsizes(EPElNino_anom(0,:))
    Comp_avg = new((/3, nlat, nlon/),"float")
    Comp_prob = new((/3, nlat, nlon/),"float")
    Comp_avg(0,:,:) = EPElNino_anom*0.9
	Comp_avg(1,:,:) = CPElNino_anom*0.9
	Comp_avg(2,:,:) = LaNina_anom*0.9
    Comp_prob(0,:,:) = EPElNino_prob*0.9
	Comp_prob(1,:,:) = CPElNino_prob*0.9;
	Comp_prob(2,:,:) = LaNina_prob*0.9

	Comp_avg!0 = "type"
	Comp_avg!1 = "lat"
	Comp_avg!2 = "lon"
	Comp_avg&lat = EPElNino_anom&lat
	Comp_avg&lon = EPElNino_anom&lon
	copy_VarCoords(Comp_avg,Comp_prob)

wks = gsn_open_wks(pltTpye, pltName)

	res = True
	res@gsnDraw = False
	res@gsnFrame = False
	res@gsnAddCyclic = True
	res@gsnMaximize = False
	res@lbLabelBarOn         = False
	res@cnFillOn = True
	res@cnLinesOn = False
	res@cnLineLabelsOn = False
	res@cnInfoLabelOn = False
	res@tmXTOn = False
	res@tmYROn = False

	;------plot prob	
	res_prob = res
	res_prob@lbLabelBarOn = False
	res_prob@cnLevelSelectionMode= "ExplicitLevels"
	res_prob@cnLevels              = (/siglvl/)
	res_prob@cnMonoFillPattern = False 
	res_prob@cnMonoFillColor = False
	res_prob@cnFillColors = (/"white","white"/)
	res_prob@cnFillPatterns        = (/-1,17/) ;
	res_prob@cnFillDotSizeF = 0.004 ; dots size
    res_prob@cnFillScaleF = 1. ;dots density
	
    res@mpMinLatF = latMin              ; 最小纬度（北极区域）
    res@mpMaxLatF = latMax              ; 最大纬度（北极）
    res@mpMinLonF = lonMin            ; 最小经度
    res@mpMaxLonF = lonMax             ; 最大经度
    res@mpCenterLonF = (lonMin+lonMax)/2
    res@mpLandFillColor = "white" ; 设置地图填充颜色
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/-12,-10,-8,-6,-4,-2,2,4,6,8,10,12/)
    res@cnFillPalette = "temp_diff_18lev"
	FontHeightF = 0.02
    res@tmYLLabelFontHeightF = FontHeightF
    res@tmXBLabelFontHeightF = FontHeightF
    res@gsnMajorLonSpacing = 60
    res@gsnLeftStringOrthogonalPosF = 0.01
	res@gsnLeftStringFontHeightF = FontHeightF
    res@gsnRightStringOrthogonalPosF = 0.01
	res@gsnRightStringFontHeightF = FontHeightF
	
	plot = new(3, graphic)
	plot_prob = new(3, graphic)
	box_line = new(3, graphic)
	LeftString = (/"(g) Cyclone day frequency (%)", "(h) Cyclone day frequency (%)", "(i) Cyclone day frequency (%)"/)
	ntilde  = "n~H-13V2F35~D~FV-2H3~"
	RightString = (/"EP El Ni"+ntilde+"o non-pattern days", "CP El Ni"+ntilde+"o non-pattern days", "La Ni"+ntilde+"a non-pattern days"/)
	do i = 0, 2
		res@gsnLeftString = LeftString(i)
		res@gsnRightString = RightString(i)
		plot(i) = gsn_csm_contour_map_ce(wks,Comp_avg(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
		
		lon1 = 153
		lon2 = 170
		lat1 = -45
		lat2 = -30
		;-----region box
		X = (/ lon1, lon2, lon2, lon1, lon1 /)
		Y = (/ lat1, lat1, lat2, lat2, lat1 /)

		pres_line = True
		pres_line@gsLineColor = "black"
		pres_line@gsLineThicknessF = 3.0
		pres_line@gsLineDashPattern = 0    ; 实线
		box_line(i) = gsn_add_polyline(wks, plot(i), X, Y, pres_line)
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 1.
	resP@pmLabelBarHeightF = 0.07
	; resP@pmLabelBarParallelPosF = 0.015
	resP@lbLabelFontHeightF = 0.02
	gsn_panel(wks,plot,(/3,1/),resP)  