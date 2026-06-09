
sphere = "SH"
latMax = 0
latMin = -90;-60
lonMax = 360;180
lonMin = 0;90
siglvl = 0.95

threshold = 6
days_window = 30+31+30

year_start = 1970
year_end = 2022
total_year = ispan(year_start, year_end, 1)
season = "SON"
pltTpye = "pdf"
pltName = "Fig7_SON"

; EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016/)-1 ;2024/)
; CP_ElNino_year = (/1987, 1992, 1995, 2003, 2005, 2010, 2019/)-1
; ; EP_LaNina_year = (/1985, 1996, 2000, 2008, 2018, 2021/)
; LaNina_year = (/1984, 1989, 1999, 2001, 2011, 2012, 2022/)-1; 2023/)

EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016/)-1;, 2024/)
CP_ElNino_year = (/1987, 1992, 1995, 2003, 2005, 2010, 2019/)-1
; LaNina_year = (/1984, 1989, 1999, 2001, 2011, 2012, 2022, 2023/)
EP_LaNina_year = (/1996, 2000, 2008, 2021/)-1
CP_LaNina_year = (/1984, 1989, 1999, 2011, 2012, 2022, 2023/)-1
LaNina_year = array_append_record(EP_LaNina_year,CP_LaNina_year,0)

f = addfile("/g/data/if69/ls3248/data_HF/weather/cyclone/cyclone."+season+".mon.anom_slidingClim.nc", "r") ;total days in each month
cyclone_temp = f->cyclone
cyclone = cyclone_temp*1. ;转换为浮点型（每月气旋总数）
copy_VarCoords(cyclone_temp, cyclone)

time = cyclone&time
lat = cyclone&lat
lon = cyclone&lon
ntime = dimsizes(time)
nlat = dimsizes(lat)
nlon = dimsizes(lon)
cyclone_season = dim_sum_n(reshape(cyclone,(/ntime/3,3,nlat,nlon/)), 1)/days_window*100.
cyclone_season!0 = "time"
cyclone_season!1 = "lat"
cyclone_season!2 = "lon"
cyclone_season&time = total_year
cyclone_season&lat = lat
cyclone_season&lon = lon

;Comp & prob
;-----Composite and T-test
	iflag    = True      
	tval_opt = False 
	;-----total
    total_avg = dim_avg_n_Wrap(cyclone_season({total_year},:,:),0)
    total_var = dim_variance_n_Wrap(cyclone_season({total_year},:,:),0)
    total_num = dimsizes(total_year)

	;-----EP_ElNino
    EP_ElNino_avg = dim_avg_n_Wrap(cyclone_season({EP_ElNino_year},:,:),0)
    EP_ElNino_var = dim_variance_n_Wrap(cyclone_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_ElNino_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(EP_ElNino_avg,EP_ElNino_prob)

	;-----CP_ElNino
    CP_ElNino_avg = dim_avg_n_Wrap(cyclone_season({CP_ElNino_year},:,:),0)
    CP_ElNino_var = dim_variance_n_Wrap(cyclone_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_ElNino_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(CP_ElNino_avg,CP_ElNino_prob)
	
	;-----LaNina
    LaNina_avg = dim_avg_n_Wrap(cyclone_season({LaNina_year},:,:),0)
    LaNina_var = dim_variance_n_Wrap(cyclone_season({LaNina_year},:,:),0)
    LaNina_num = dimsizes(LaNina_year)
    LaNina_prob = 1. - ttest(LaNina_avg, LaNina_var, LaNina_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(LaNina_avg,LaNina_prob)

    ;-----comp
    Comp_avg = new((/3, nlat, nlon/),"float")
    Comp_prob = new((/3, nlat, nlon/),"float")
    Comp_avg(0,:,:) = EP_ElNino_avg
	Comp_avg(1,:,:) = CP_ElNino_avg;
	Comp_avg(2,:,:) = LaNina_avg
    Comp_prob(0,:,:) = EP_ElNino_prob
	Comp_prob(1,:,:) = CP_ElNino_prob;
	Comp_prob(2,:,:) = LaNina_prob

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
	LeftString = (/"(a) Cyclone day frequency (%)", "(b) Cyclone day frequency (%)", "(c) Cyclone day frequency (%)"/)
	ntilde  = "n~H-13V2F35~D~FV-2H3~"
	RightString = (/"EP El Ni"+ntilde+"o SON", "CP El Ni"+ntilde+"o SON", "La Ni"+ntilde+"a SON"/)
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