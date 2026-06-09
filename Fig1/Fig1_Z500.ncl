year_start = 1980
year_end = 2024
total_year = ispan(year_start,year_end,1)
level = 500

latS = -90
latN = 0
lonE = 0
lonW = 360

season = "SON"
FontHeightF = 0.025

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
    
siglvl = 0.95

pltType = "pdf"
pltName = "Fig1_Z"+level

;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/hgt/hgt.ERA5.mon."+level+"hPa.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	hgt = lonFlip(rm_single_dims(f->hgt(start_ind:end_ind,:,:,:)))
	hgt_season = rm_single_dims(month_to_seasonN(hgt,season))
    nlat = dimsizes(hgt&latitude)
    nlon = dimsizes(hgt&longitude)
    hgt_season&time = total_year
    
;-----Composite and T-test
	iflag    = True      
	tval_opt = False 
	;-----total
    total_avg = dim_avg_n_Wrap(hgt_season({total_year},:,:),0)
    total_var = dim_variance_n_Wrap(hgt_season({total_year},:,:),0)
    total_num = dimsizes(total_year)

	;-----EP_ElNino
    EP_ElNino_avg = dim_avg_n_Wrap(hgt_season({EP_ElNino_year},:,:),0)
    EP_ElNino_var = dim_variance_n_Wrap(hgt_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_ElNino_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(EP_ElNino_avg,EP_ElNino_prob)

	;-----CP_ElNino
    CP_ElNino_avg = dim_avg_n_Wrap(hgt_season({CP_ElNino_year},:,:),0)
    CP_ElNino_var = dim_variance_n_Wrap(hgt_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_ElNino_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(CP_ElNino_avg,CP_ElNino_prob)
	
	;-----LaNina
    LaNina_avg = dim_avg_n_Wrap(hgt_season({LaNina_year},:,:),0)
    LaNina_var = dim_variance_n_Wrap(hgt_season({LaNina_year},:,:),0)
    LaNina_num = dimsizes(LaNina_year)
    LaNina_prob = 1. - ttest(LaNina_avg, LaNina_var, LaNina_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(LaNina_avg,LaNina_prob)

    ;-----combine
    Comp_avg = new((/3, nlat, nlon/),"float")
    Comp_prob = new((/3, nlat, nlon/),"float")
    Comp_avg(0,:,:) = EP_ElNino_avg
	Comp_avg(1,:,:) = CP_ElNino_avg;
	Comp_avg(2,:,:) = LaNina_avg
    Comp_prob(0,:,:) = EP_ElNino_prob
	Comp_prob(1,:,:) = CP_ElNino_prob;
	Comp_prob(2,:,:) = LaNina_prob
    
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
	
	res@mpCenterLonF = (lonE+lonW)/2.
    res@mpMinLatF = latS
    res@mpMaxLatF = latN
	res@mpMinLonF = lonE
	res@mpMaxLonF = lonW
	res@mpShapeMode = "FreeAspect"
	; res@vpHeightF         = 0.5                    ; change aspect ratio of plot
	; res@vpWidthF          = 1.
	; res@cnLevelSelectionMode = "ManualLevels"
	; res@cnLevelSpacingF = 5
	; res@cnMinLevelValF = -40
	; res@cnMaxLevelValF = 40
	; res@cnFillPalette = "BlueDarkRed18"
	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLevels = (/-50,-40,-30,-20,-15,-10,-5,5,10,15,20,30,40,50/)
	res@cnFillPalette = "NCV_blu_red"
	res@tmYLLabelFontHeightF = FontHeightF
	res@tmXBLabelFontHeightF = FontHeightF
	res@gsnMajorLonSpacing = 60

	res@gsnLeftStringOrthogonalPosF = 0.01
	res@gsnRightStringOrthogonalPosF = 0.01
	res@gsnRightStringFontHeightF = FontHeightF
	res@gsnLeftStringFontHeightF = FontHeightF

	ntilde  = "n~H-13V2F35~D~FV-2H3~"
	RightString = (/"EP El Ni"+ntilde+"o","CP El Ni"+ntilde+"o","La Ni"+ntilde+"a"/)
	LeftString = (/"(g) Z500","(h) Z500","(i) Z500"/)
	
	plot = new(3, graphic)
	plot_prob = new(3, graphic)
	ploy1 = new(3, graphic)
	ploy2 = new(3, graphic)
	do i = 0, 2
		res@gsnLeftString = LeftString(i)
		res@gsnRightString = RightString(i)
		plot(i) = gsn_csm_contour_map_ce(wks,Comp_avg(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))

		; Add the rectangular box
		box_res = True
		box_res@gsLineColor = "red"
		box_res@gsLineThicknessF = 2.0
		box_res@gsLineDashPattern = 0

		; Define box coordinates
		latS_box = -80.0  ; 80S
		latN_box = -20.0  ; 20S
		lonW_box = 280.  ; 80W
		lonE_box = 80.0  ;80E
		lonC = 180.
		lon_box1 = (/lonC, lonE_box, lonE_box, lonC/)
		lat_box1 = (/latN_box, latN_box, latS_box, latS_box/)

		lon_box2 = (/lonC, lonW_box, lonW_box, lonC/)
		lat_box2 = (/latN_box, latN_box, latS_box, latS_box/)
		; Draw the box
		ploy1(i) = gsn_add_polyline(wks, plot(i), lon_box1, lat_box1, box_res)
		ploy2(i) = gsn_add_polyline(wks, plot(i), lon_box2, lat_box2, box_res)
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.65
	resP@pmLabelBarHeightF = 0.04
	resP@lbLabelFontHeightF = 0.01
	gsn_panel(wks,plot,(/1,3/),resP)  

	; ;-----output the Comp_avg
	; out_file = "composite_pattern.nc"
	; if (fileexists(out_file)) then
	; 	system("rm -f "+out_file)
	; end if
	; f_out = addfile(out_file,"c")
	; Comp_avg!0 = "ENSO_type"
	; Comp_avg&ENSO_type = (/0,1,2/) ;0-EPElNino 1-CPElNino 2-LaNina
	; f_out->composite_pattern = Comp_avg

	; print(pattern_cor(Comp_avg(0,:,:), Comp_avg(1,:,:), 1., 0))
	; print(pattern_cor(Comp_avg(0,:,:), Comp_avg(2,:,:), 1., 0))
	; print(pattern_cor(Comp_avg(1,:,:), Comp_avg(2,:,:), 1., 0))