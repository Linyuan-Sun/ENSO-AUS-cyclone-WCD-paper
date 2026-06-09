type = "pattern" ;"pattern" "nonpattern" "season"
pltType = "pdf"
pltName = "Fig12_"+type

siglvl = (/0.95, 0.95/)
FontHeightF = 0.025

latMin = -44.5
latMax = -10
lonMin = 112
lonMax = 156

var_EP = "var_EPElNino_"+type+"_mean"
var_CP = "var_CPElNino_"+type+"_mean"
var_La = "var_LaNina_"+type+"_mean"

f = addfile("rainfall_season_pattern_nonpattern.nc","r")
EP = f->$var_EP$({latMin:latMax},{lonMin:lonMax})
CP = f->$var_CP$({latMin:latMax},{lonMin:lonMax})
La = f->$var_La$({latMin:latMax},{lonMin:lonMax})

lat = EP&lat
lon = EP&lon
nlat = dimsizes(lat)
nlon = dimsizes(lon)

Comp_mean = new((/3,nlat,nlon/),"float")
Comp_mean(0,:,:) = EP*1.3
Comp_mean(1,:,:) = CP*1.1
Comp_mean(2,:,:) = La*1.1
Comp_mean!0 = "type"
Comp_mean!1 = "lat"
Comp_mean!2 = "lon"
Comp_mean&lat = lat
Comp_mean&lon = lon

;prob
var_EP = "var_EPElNino_"+type+"_prob"
var_CP = "var_CPElNino_"+type+"_prob"
var_La = "var_LaNina_"+type+"_prob"
Comp_prob = new((/3, nlat, nlon/),"float")
Comp_prob(0,:,:) = f->$var_EP$({latMin:latMax},{lonMin:lonMax})
Comp_prob(1,:,:) = f->$var_CP$({latMin:latMax},{lonMin:lonMax})
Comp_prob(2,:,:) = f->$var_La$({latMin:latMax},{lonMin:lonMax})
copy_VarCoords(Comp_mean,Comp_prob)

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
	res_prob@cnMonoFillPattern = False 
	res_prob@cnMonoFillColor = False
	res_prob@cnFillColors = (/"white","white"/)
	res_prob@cnFillPatterns        = (/-1,17/) ;
	res_prob@cnFillDotSizeF = 0.004 ; dots size
    res_prob@cnFillScaleF = 1.2

	res@mpMinLatF = latMin 
    res@mpMaxLatF = latMax
    res@mpMinLonF = lonMin   
    res@mpMaxLonF = lonMax
    res@mpCenterLonF = (lonMin+lonMax)/2.
	res@mpShapeMode = "FreeAspect"
	; res@vpHeightF         = 0.5                    ; change aspect ratio of plot
	; res@vpWidthF          = 1.5
	res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1/);
	cmap = read_colormap_file("CBR_drywet")
	cmap(5,:) = (/0,0,0,0/)
    res@cnFillPalette = cmap
	res@tmYLLabelFontHeightF = FontHeightF
	res@tmXBLabelFontHeightF = FontHeightF
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
	res@gsnMajorLonSpacing = 15

	res@gsnLeftStringOrthogonalPosF = 0.01
	res@gsnRightStringOrthogonalPosF = 0.01
	res@gsnRightStringFontHeightF = FontHeightF
	res@gsnLeftStringFontHeightF = FontHeightF
	ntilde  = "n~H-13V2F35~D~FV-2H3~"
	rightString = (/"EP El Ni"+ntilde+"o pattern","CP El Ni"+ntilde+"o pattern","La Ni"+ntilde+"a pattern"/)
	leftString = (/"(d) rainfall","(e) rainfall","(f) rainfall"/)
	
	plot = new(3, graphic)
	plot_prob = new(3, graphic)
    do i = 0, 2
		res@gsnLeftString = leftString(i)
		res@gsnRightString = rightString(i)
		if i.eq.0 then
			res_prob@cnLevels              = (/siglvl(0)/)
		else
			res_prob@cnLevels              = (/siglvl(1)/)
		end if
		plot(i) = gsn_csm_contour_map_ce(wks,Comp_mean(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.6
	resP@pmLabelBarHeightF = 0.045
	; resP@pmLabelBarParallelPosF = 0.015
	resP@lbLabelFontHeightF = 0.015
	gsn_panel(wks,plot,(/3,1/),resP)  


