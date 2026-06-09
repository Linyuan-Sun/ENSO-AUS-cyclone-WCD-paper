pattern = "z500"
season = "SON"
var = "z"
latS = -90.
latN = 0.
lonW = 40.
lonE = 320.
level = 500

delta_day = 10
plot_day_num = (/0,1/)
plot_day = plot_day_num + delta_day
nplot_day = dimsizes(plot_day)

type = (/"EP ElNino","CP ElNino","LaNina"/)
ntype = dimsizes(type)

pltType = "pdf"
pltName = pattern+"_pattern_"+var+level+"_plot"

siglvl = 0.95

; load data
f = addfile(pattern+"_pattern_"+var+level+"_"+season+".nc", "r")
varName = var+"_EPElNino"
var1 = f->$varName$(:,plot_day,{latS:latN},{lonW:lonE})
varName = var+"_CPElNino"
var2 = f->$varName$(:,plot_day,{latS:latN},{lonW:lonE})
varName = var+"_LaNina"
var3 = f->$varName$(:,plot_day,{latS:latN},{lonW:lonE})

; Comp and T-test
iflag    = True      
tval_opt = False 

var1_avg = dim_avg_n_Wrap(var1,0)
var1_std = dim_stddev_n_Wrap(var1,0)
n = dimsizes(var1(:,0,0,0))
var1_tvalue = var1_avg/(var1_std/sqrt(n))
df = n - 1
var1_prob = 1 - student_t(var1_tvalue, df)
copy_VarCoords(var1_avg,var1_prob)

var2_avg = dim_avg_n_Wrap(var2,0)
var2_std = dim_stddev_n_Wrap(var2,0)
n = dimsizes(var2(:,0,0,0))
var2_tvalue = var2_avg/(var2_std/sqrt(n))
df = n - 1
var2_prob = 1 - student_t(var2_tvalue, df)
copy_VarCoords(var2_avg,var2_prob)

var3_avg = dim_avg_n_Wrap(var3,0)
var3_std = dim_stddev_n_Wrap(var3,0)
n = dimsizes(var3(:,0,0,0))
var3_tvalue = var3_avg/(var3_std/sqrt(n))
df = n - 1
var3_prob = 1 - student_t(var3_tvalue, df)
copy_VarCoords(var3_avg,var3_prob)

; var_plot
if var.eq."z" then
    var1_plot = var1_avg/9.8
    var2_plot = var2_avg/9.8
    var3_plot = var3_avg/9.8
else
    var1_plot = var1_avg
    var2_plot = var2_avg
    var3_plot = var3_avg
end if

lat = var1_avg&lat
lon = var1_avg&lon
nday = nplot_day
nlat = dimsizes(lat)
nlon = dimsizes(lon)
var_plot = new((/ntype,nday,nlat,nlon/),"float")
var_plot(0,:,:,:) = var1_plot
var_plot(1,:,:,:) = var2_plot
var_plot(2,:,:,:) = var3_plot
var_plot!0 = "type"
var_plot!1 = "day"
var_plot!2 = "lat"
var_plot!3 = "lon"
var_plot&lat = lat
var_plot&lon = lon

var_prob = new((/ntype,nday,nlat,nlon/),"float")
var_prob(0,:,:,:) = var1_prob
var_prob(1,:,:,:) = var2_prob
var_prob(2,:,:,:) = var3_prob
copy_VarCoords(var_plot,var_prob)

var_plot_mean = dim_avg_n_Wrap(var_plot,1)
var_prob_mean = dim_avg_n_Wrap(var_prob,1)

;-----plot figure
	wks = gsn_open_wks(pltType,pltName)
	res = True
	res@gsnDraw = False
	res@gsnFrame = False
	res@gsnAddCyclic = False
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
    res_prob@cnFillScaleF = 0.6 ;dots density
	
	res@mpCenterLonF = (lonW+lonE)/2.
	res@mpMinLatF = -90
	res@mpMaxLatF = 0
	res@mpMinLonF = 40
	res@mpMaxLonF = 320
	; res@mpShapeMode = "FreeAspect"
	; res@vpHeightF         = 0.5                    ; change aspect ratio of plot
	; res@vpWidthF          = 1.
	; res@cnLevelSelectionMode = "ManualLevels"
	; res@cnLevelSpacingF = 5
	; res@cnMinLevelValF = -40
	; res@cnMaxLevelValF = 40
	; res@cnFillPalette = "BlueDarkRed18"
	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLevels = (/-50,-40,-30,-20,-10,-5,5,10,20,30,40,50/)
	res@cnFillPalette = "NCV_blu_red"
	res@tmYLLabelFontHeightF = 0.025
	res@tmXBLabelFontHeightF = 0.025
	res@gsnMajorLonSpacing = 60
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnRightString = "";"gpm"
	res@gsnRightStringOrthogonalPosF = 0.
	
	plot = new(3, graphic)
	plot_prob = new(3, graphic)
	ploy1 = new(3, graphic)
	ploy2 = new(3, graphic)
	; LeftString = (/"(a) EP ElNino "+season+": Z"+level,"(b) CP ElNino "+season+": Z"+level,\
	; 			   "(c) LaNina "+season+": Z"+level/)
	do i = 0, 2
		res@gsnLeftString = "";LeftString(i)
		plot(i) = gsn_csm_contour_map_ce(wks,var_plot_mean(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,var_prob_mean(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))

		; Add the rectangular box
		box_res = True
		box_res@gsLineColor = "red"
		box_res@gsLineThicknessF = 2.0
		box_res@gsLineDashPattern = 0

		; Define box coordinates
		latS = -80.0  ; 80S
		latN = -20.0  ; 20S
		lonW = 280.  ; 80W
		lonE = 80.0  ;80E
		lonC = 180.
		lon_box1 = (/lonC, lonE, lonE, lonC/)
		lat_box1 = (/latN, latN, latS, latS/)

		lon_box2 = (/lonC, lonW, lonW, lonC/)
		lat_box2 = (/latN, latN, latS, latS/)
		; Draw the box
		ploy1(i) = gsn_add_polyline(wks, plot(i), lon_box1, lat_box1, box_res)
		ploy2(i) = gsn_add_polyline(wks, plot(i), lon_box2, lat_box2, box_res)
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.5
	resP@pmLabelBarHeightF = 0.045
	resP@lbLabelFontHeightF = 0.01
	gsn_panel(wks,plot,(/1,3/),resP)  