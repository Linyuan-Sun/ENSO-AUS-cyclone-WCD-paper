year_start = 1980
year_end = 2024
total_year = ispan(year_start,year_end,1)

latS = -30
latN = 30
lonE = 30
lonW = 270

season = "SON"

;-----SON define
; EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016, 2024/)
; CP_ElNino_year = (/1987, 1992, 1995, 2003, 2005, 2010, 2019/)
; EP_LaNina_year = (/1985, 1996, 2000, 2008, 2018, 2021/)
; CP_LaNina_year = (/1984, 1989, 1999, 2001, 2011, 2012, 2022, 2023/)
; LaNina_year = array_append_record(EP_LaNina_year,CP_LaNina_year,0)
EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016, 2024/)
CP_ElNino_year = (/1987, 1992, 1995, 2003, 2005, 2010, 2019/)
EP_LaNina_year = (/1996, 2000, 2008, 2021/)
CP_LaNina_year = (/1984, 1989, 1999, 2011, 2012, 2022, 2023/)
LaNina_year = array_append_record(EP_LaNina_year,CP_LaNina_year,0)

pltType = "pdf"
pltName = "Fig1_OLR"
FontHeightF = 0.025

siglvl = 0.95

;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/OLR/olr.mon.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	olr = lonFlip(f->olr(start_ind:end_ind,{latS:latN},:)) ;make lon consistent with SST
	olr_season = rm_single_dims(month_to_seasonN(olr,season))
	olr_season&time = total_year
    nlat_olr = dimsizes(olr&lat)
    nlon_olr = dimsizes(olr&lon)
	delete(time)

;-----Composite and T-test
	iflag    = True      
	tval_opt = False 

	;-----olr
	;-----total
	total_avg_olr = dim_avg_n_Wrap(olr_season({total_year},:,:),0)
	total_var_olr = dim_variance_n_Wrap(olr_season({total_year},:,:),0)
	total_num_olr = dimsizes(total_year)
	;-----EP ElNino
    EP_ElNino_avg_olr = dim_avg_n_Wrap(olr_season({EP_ElNino_year},:,:),0)
    EP_ElNino_var_olr = dim_variance_n_Wrap(olr_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num_olr = dimsizes(EP_ElNino_year)
    EP_ElNino_prob_olr = 1. - ttest(EP_ElNino_avg_olr, EP_ElNino_var_olr, EP_ElNino_num_olr, total_avg_olr, total_var_olr, total_num_olr, iflag, tval_opt)
	copy_VarCoords(EP_ElNino_avg_olr,EP_ElNino_prob_olr)
	;-----CP ElNino
    CP_ElNino_avg_olr = dim_avg_n_Wrap(olr_season({CP_ElNino_year},:,:),0)
    CP_ElNino_var_olr = dim_variance_n_Wrap(olr_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num_olr = dimsizes(CP_ElNino_year)
    CP_ElNino_prob_olr = 1. - ttest(CP_ElNino_avg_olr, CP_ElNino_var_olr, CP_ElNino_num_olr, total_avg_olr, total_var_olr, total_num_olr, iflag, tval_opt)
	copy_VarCoords(CP_ElNino_avg_olr,CP_ElNino_prob_olr)
	;-----LaNina
	LaNina_avg_olr = dim_avg_n_Wrap(olr_season({LaNina_year},:,:),0)
	LaNina_var_olr = dim_variance_n_Wrap(olr_season({LaNina_year},:,:),0)
	LaNina_num_olr = dimsizes(LaNina_year)
	LaNina_prob_olr = 1. - ttest(LaNina_avg_olr, LaNina_var_olr, LaNina_num_olr, total_avg_olr, total_var_olr, total_num_olr, iflag, tval_opt)
	copy_VarCoords(LaNina_avg_olr,LaNina_prob_olr)
    ;-----combine
    Comp_avg_olr = new((/3, nlat_olr, nlon_olr/),"float")
    Comp_prob_olr = new((/3, nlat_olr, nlon_olr/),"float")
    Comp_avg_olr(0,:,:) = EP_ElNino_avg_olr
    Comp_avg_olr(1,:,:) = CP_ElNino_avg_olr
	Comp_avg_olr(2,:,:) = LaNina_avg_olr
    Comp_prob_olr(0,:,:) = EP_ElNino_prob_olr
    Comp_prob_olr(1,:,:) = CP_ElNino_prob_olr
	Comp_prob_olr(2,:,:) = LaNina_prob_olr

;-----plot figure
	wks = gsn_open_wks(pltType,pltName)
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
	res_prob@cnFillDotSizeF = 0.0025 ; dots size
    res_prob@cnFillScaleF = 0.65 ;dots density
	; res_prob@cnFillDrawOrder = "PostDraw" ; dots on top of contour

	res@mpShapeMode = "FreeAspect"
	res@mpCenterLonF = (lonE+lonW)/2.
    res@mpMinLatF = latS
    res@mpMaxLatF = latN
	res@mpMinLonF = lonE
	res@mpMaxLonF = lonW
	cmap = read_colormap_file("MPL_BrBG")
	res@mpLandFillColor = cmap(65-2,:)
	res@mpGeophysicalLineColor = "black"
	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLevels = (/-30,-25,-20,-15,-10,-5,5,10,15,20,25,30/)
	res@cnFillPalette = "NCV_blu_red"
	res@tmYLLabelFontHeightF = FontHeightF
	res@tmXBLabelFontHeightF = FontHeightF
	res@gsnMajorLonSpacing = 60
	res@gsnMajorLatSpacing = 20
	res@gsnLeftStringOrthogonalPosF = 0.01
	res@gsnRightStringOrthogonalPosF = 0.01
	res@gsnRightStringFontHeightF = FontHeightF
	res@gsnLeftStringFontHeightF = FontHeightF
	res@cnSmoothingOn = True
	
	ntilde  = "n~H-13V2F35~D~FV-2H3~"
	RightString = (/"EP El Ni"+ntilde+"o","CP El Ni"+ntilde+"o","La Ni"+ntilde+"a"/)
	LeftString = (/"(d) OLR","(e) OLR","(f) OLR"/)

	plot = new(3, graphic)
	plot_prob = new(3, graphic)
	ploy1 = new(3, graphic)
	ploy2 = new(3, graphic)
	do i = 0,2

		res@gsnLeftString = LeftString(i)
		res@gsnRightString = RightString(i)
		plot(i) = gsn_csm_contour_map_ce(wks,Comp_avg_olr(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob_olr(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))

		; ; Add the rectangular box
		; box_res = True
		; box_res@gsLineColor = "red"
		; box_res@gsLineThicknessF = 2.0
		; box_res@gsLineDashPattern = 0

		; ; Define box coordinates
		; latS_box = -10.0  ; 10S
		; latN_box = 10.0  ; 10N
		; lonW_box = 210.  ; 150W
		; lonE_box = 80.0  ;80E
		; lonC = 180.
		; lon_box1 = (/lonC, lonE_box, lonE_box, lonC/)
		; lat_box1 = (/latN_box, latN_box, latS_box, latS_box/)

		; lon_box2 = (/lonC, lonW_box, lonW_box, lonC/)
		; lat_box2 = (/latN_box, latN_box, latS_box, latS_box/)
		; ; Draw the box
		; ploy1(i) = gsn_add_polyline(wks, plot(i), lon_box1, lat_box1, box_res)
		; ploy2(i) = gsn_add_polyline(wks, plot(i), lon_box2, lat_box2, box_res)
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.65
	resP@pmLabelBarHeightF = 0.04
	resP@lbLabelFontHeightF = 0.01
	gsn_panel(wks,plot,(/1,3/),resP)  