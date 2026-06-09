year_start = 1980
year_end = 2024
total_year = ispan(year_start,year_end,1)

latS = -35
latN = 35
lonE = 30
lonW = 290

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
pltName = "Fig1_SST"
FontHeightF = 0.025

siglvl = 0.95

;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/HadISST/HadISST.mon.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	sst = f->sst(start_ind:end_ind,{latS:latN},:)
	sst_season = rm_single_dims(month_to_seasonN(sst,season))
	sst_season&time = total_year
    nlat_sst = dimsizes(sst&latitude)
    nlon_sst = dimsizes(sst&longitude)
	delete(time)
    
;-----Composite and T-test
	iflag    = True      
	tval_opt = False 
	
	;-----sst
	;-----total
	total_avg_sst = dim_avg_n_Wrap(sst_season({total_year},:,:),0)
	total_var_sst = dim_variance_n_Wrap(sst_season({total_year},:,:),0)
	total_num_sst = dimsizes(total_year)
	;-----EP ElNino
    EP_ElNino_avg_sst = dim_avg_n_Wrap(sst_season({EP_ElNino_year},:,:),0)
    EP_ElNino_var_sst = dim_variance_n_Wrap(sst_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num_sst = dimsizes(EP_ElNino_year)
    EP_ElNino_prob_sst = 1. - ttest(EP_ElNino_avg_sst, EP_ElNino_var_sst, EP_ElNino_num_sst, total_avg_sst, total_var_sst, total_num_sst, iflag, tval_opt)
	copy_VarCoords(EP_ElNino_avg_sst,EP_ElNino_prob_sst)
	;-----CP ElNino
    CP_ElNino_avg_sst = dim_avg_n_Wrap(sst_season({CP_ElNino_year},:,:),0)
    CP_ElNino_var_sst = dim_variance_n_Wrap(sst_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num_sst = dimsizes(CP_ElNino_year)
    CP_ElNino_prob_sst = 1. - ttest(CP_ElNino_avg_sst, CP_ElNino_var_sst, CP_ElNino_num_sst, total_avg_sst, total_var_sst, total_num_sst, iflag, tval_opt)
	copy_VarCoords(CP_ElNino_avg_sst,CP_ElNino_prob_sst)
	;-----LaNina
	LaNina_avg_sst = dim_avg_n_Wrap(sst_season({LaNina_year},:,:),0)
	LaNina_var_sst = dim_variance_n_Wrap(sst_season({LaNina_year},:,:),0)
	LaNina_num_sst = dimsizes(LaNina_year)
	LaNina_prob_sst = 1. - ttest(LaNina_avg_sst, LaNina_var_sst, LaNina_num_sst, total_avg_sst, total_var_sst, total_num_sst, iflag, tval_opt)
	copy_VarCoords(LaNina_avg_sst,LaNina_prob_sst)
    ;-----combine
    Comp_avg_sst = new((/3, nlat_sst, nlon_sst/),"float")
    Comp_prob_sst = new((/3, nlat_sst, nlon_sst/),"float")
    Comp_avg_sst(0,:,:) = EP_ElNino_avg_sst
    Comp_avg_sst(1,:,:) = CP_ElNino_avg_sst
	Comp_avg_sst(2,:,:) = LaNina_avg_sst
    Comp_prob_sst(0,:,:) = EP_ElNino_prob_sst
    Comp_prob_sst(1,:,:) = CP_ElNino_prob_sst
	Comp_prob_sst(2,:,:) = LaNina_prob_sst

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

    res@mpMinLatF = latS
    res@mpMaxLatF = latN
	res@mpMinLonF = lonE
	res@mpMaxLonF = lonW
	res@mpCenterLonF = (lonE+lonW)/2.
	cmap = read_colormap_file("MPL_BrBG")
	res@mpLandFillColor = cmap(65-2,:)
	res@mpGeophysicalLineColor = "grey50"
	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLevels = (/-1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5/)
	res@cnFillPalette = "NCV_blu_red"
	res@tmYLLabelFontHeightF = FontHeightF
	res@tmXBLabelFontHeightF = FontHeightF
	res@gsnMajorLonSpacing = 60
	res@gsnMajorLatSpacing = 20
	res@gsnLeftStringOrthogonalPosF = 0.01
	res@gsnRightStringOrthogonalPosF = 0.01
	res@gsnRightStringFontHeightF = FontHeightF
	res@gsnLeftStringFontHeightF = FontHeightF
	
	ntilde  = "n~H-13V2F35~D~FV-2H3~"
	RightString = (/"EP El Ni"+ntilde+"o","CP El Ni"+ntilde+"o","La Ni"+ntilde+"a"/)
	LeftString = (/"(a) SST","(b) SST","(c) SST"/)

	plot = new(3, graphic)
	plot_prob = new(3, graphic)
	do i = 0,2
		res@gsnLeftString = LeftString(i)
		res@gsnRightString = RightString(i)
		plot(i) = gsn_csm_contour_map_ce(wks,Comp_avg_sst(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob_sst(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.65
	resP@pmLabelBarHeightF = 0.04
	resP@lbLabelFontHeightF = 0.01
	gsn_panel(wks,plot,(/1,3/),resP)  