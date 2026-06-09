type = "seasonal" 
pltType = "pdf"
pltName = "Fig3_"+type

siglvl = 0.95

latMin = -80
latMax = -20;25
lonMin = 80
lonMax = 280

FontHeightF = 0.045

year_start = 1980
year_end = 2024
total_year = ispan(year_start,year_end,1)
EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016, 2024/)
CP_ElNino_year = (/1987, 1992, 1995, 2003, 2005, 2010, 2019/)
EP_LaNina_year = (/1996, 2000, 2008, 2021/)
CP_LaNina_year = (/1984, 1989, 1999, 2011, 2012, 2022, 2023/)
LaNina_year = array_append_record(EP_LaNina_year,CP_LaNina_year,0)

season = "SON"
level = 500

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
	res@gsnAddCyclic = False
	res@gsnMaximize = False
	res@cnFillOn = True
	res@cnLinesOn = False
	res@cnLineLabelsOn = False
	res@cnInfoLabelOn = False
    res@lbLabelBarOn         = False
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
	res_prob@cnFillDotSizeF = 0.005 ; dots size
    res_prob@cnFillScaleF = 1.5

	res@mpMinLatF = latMin 
    res@mpMaxLatF = latMax
    res@mpMinLonF = lonMin   
    res@mpMaxLonF = lonMax
    res@mpCenterLonF = (lonMin+lonMax)/2.
	res@mpShapeMode = "FreeAspect"
	res@vpHeightF         = 0.5                    ; change aspect ratio of plot
	res@vpWidthF          = 1.5
	res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/-100,-80,-60,-40,-30,-20,-10,10,20,30,40,60,80,100/)
    res@cnFillPalette = "NCV_blu_red"
	res@tmYLLabelFontHeightF = FontHeightF
	res@tmXBLabelFontHeightF = FontHeightF
	res@gsnMajorLonSpacing = 60
    res@pmLabelBarWidthF = 1.5

	res@gsnLeftStringOrthogonalPosF = 0.01
	res@gsnRightStringOrthogonalPosF = 0.01
	res@gsnRightStringFontHeightF = FontHeightF
	res@gsnLeftStringFontHeightF = FontHeightF
	ntilde  = "n~H-13V2F35~D~FV-2H3~"

	rightString = (/"EP El Ni"+ntilde+"o SON (605d)","CP El Ni"+ntilde+"o SON (726d)","La Ni"+ntilde+"a SON (1331d)"/)
	leftString = (/"(a) Z500","(b) Z500","(c) Z500"/)
	
	plot = new(3, graphic)
	plot_prob = new(3, graphic)
    do i = 0, 2
		; if i.eq.2 then
		; 	res@tmXBLabelsOn = True
		; else
		; 	res@tmXBLabelsOn = False
		; end if

		res@gsnLeftString = leftString(i)
		res@gsnRightString = rightString(i)
		plot(i) = gsn_csm_contour_map_ce(wks,Comp_avg(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 1.
	resP@pmLabelBarHeightF = 0.065
	resP@lbLabelFontHeightF = 0.02
	gsn_panel(wks,plot,(/3,1/),resP)  


