
sphere = "SH"
latMax = 0
latMin = -90;-60
lonMax = 360;180
lonMin = 0;90
siglvl = 0.99

season = "SON"
pltTpye = "pdf"
pltName = "Fig8_SON"

f = addfile("season_pattern_nonpattern_out.nc", "r")
EP_mean = f->var_EPElNino_season_mean
EP_prob = f->var_EPElNino_season_prob
CP_mean = f->var_CPElNino_season_mean
CP_prob = f->var_CPElNino_season_prob
LN_mean = f->var_LaNina_season_mean
LN_prob = f->var_LaNina_season_prob

nlat = dimsizes(EP_mean(:,0))
nlon = dimsizes(EP_mean(0,:))

    ;-----comp
    Comp_avg = new((/3, nlat, nlon/),"float")
    Comp_prob = new((/3, nlat, nlon/),"float")
    Comp_avg(0,:,:) = EP_mean
	Comp_avg(1,:,:) = CP_mean;
	Comp_avg(2,:,:) = LN_mean
    Comp_prob(0,:,:) = EP_prob
	Comp_prob(1,:,:) = CP_prob
	Comp_prob(2,:,:) = LN_prob

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
    res@cnLevels = (/-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4,0.5/)
    res@cnFillPalette = "NCV_blu_red"

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
	LeftString = (/"(a) Baroclinicity", "(b) Baroclinicity", "(c) Baroclinicity"/)
	ntilde  = "n~H-13V2F35~D~FV-2H3~"
	RightString = (/"EP El Ni"+ntilde+"o SON", "CP El Ni"+ntilde+"o SON", "La Ni"+ntilde+"a SON"/)
	do i = 0, 2
		res@gsnLeftString = LeftString(i)
		res@gsnRightString = RightString(i)
		plot(i) = gsn_csm_contour_map_ce(wks,Comp_avg(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 1.
	resP@pmLabelBarHeightF = 0.07
	; resP@pmLabelBarParallelPosF = 0.015
	resP@lbLabelFontHeightF = 0.02
	gsn_panel(wks,plot,(/3,1/),resP)  