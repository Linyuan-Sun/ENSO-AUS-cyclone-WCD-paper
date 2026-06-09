type = "pattern" ;"pattern" "nonpattern" "season"
pltType = "pdf"
pltName = "Fig11_w"

siglvl = 0.95

FontHeightF = 0.03

latMin = -60
latMax = 0
lonMin = 100
lonMax = 190

var_EP = "var_EPElNino_"+type+"_mean"
var_CP = "var_CPElNino_"+type+"_mean"
var_La = "var_LaNina_"+type+"_mean"

var_EP_prob = "var_EPElNino_"+type+"_prob"
var_CP_prob = "var_CPElNino_"+type+"_prob"
var_La_prob = "var_LaNina_"+type+"_prob"

f = addfile("pattern_out_w.nc","r")
EP_mean = f->$var_EP$({latMin:latMax},{lonMin:lonMax})
CP_mean = f->$var_CP$({latMin:latMax},{lonMin:lonMax})
La_mean = f->$var_La$({latMin:latMax},{lonMin:lonMax})

EP_prob = f->$var_EP_prob$({latMin:latMax},{lonMin:lonMax})
CP_prob = f->$var_CP_prob$({latMin:latMax},{lonMin:lonMax})
La_prob = f->$var_La_prob$({latMin:latMax},{lonMin:lonMax})

lat = EP_mean&lat
lon = EP_mean&lon
nlat = dimsizes(lat)
nlon = dimsizes(lon)

Comp_mean = new((/3,nlat,nlon/),"float")
Comp_mean(0,:,:) = EP_mean*1e2
Comp_mean(1,:,:) = CP_mean*1e2
Comp_mean(2,:,:) = La_mean*1e2
Comp_mean!0 = "type"
Comp_mean!1 = "lat"
Comp_mean!2 = "lon"
Comp_mean&lat = lat
Comp_mean&lon = lon

Comp_prob = new((/3,nlat,nlon/),"float")
Comp_prob(0,:,:) = EP_prob
Comp_prob(1,:,:) = CP_prob
Comp_prob(2,:,:) = La_prob

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
	; res@vpHeightF         = 0.5                    ; change aspect ratio of plot
	; res@vpWidthF          = 1.5
	res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3/)
	cmap = read_colormap_file("MPL_RdBu")
    res@cnFillPalette = cmap(::-1,:)
	res@tmYLLabelFontHeightF = FontHeightF
	res@tmXBLabelFontHeightF = FontHeightF
	res@gsnMajorLatSpacing = 15
    res@pmLabelBarWidthF = 1.5
	res@tmXBMode = "Explicit"
	res@tmXBValues = (/100,115,130,145, 160, 175, 190/)
	res@tmXBLabels = (/"100E","115E","130E","145E","160E","175E","170W"/)
	res@tmXBMinorValues = ispan(100,190,5)
	; res@mpLandFillColor = "white"
	; res@cnFillDrawOrder      = "Predraw" 

	res@gsnLeftStringOrthogonalPosF = 0.01
	res@gsnRightStringOrthogonalPosF = 0.01
	res@gsnRightStringFontHeightF = FontHeightF
	res@gsnLeftStringFontHeightF = FontHeightF
	ntilde  = "n~H-13V2F35~D~FV-2H3~"
	rightString = (/"EP El Ni"+ntilde+"o pattern","CP El Ni"+ntilde+"o pattern","La Ni"+ntilde+"a pattern"/)
	leftString = (/"(a) w500","(b) w500","(c) w500"/)

	plot = new(3, graphic)
	plot_prob = new(3, graphic)
    do i = 0, 2
		res@gsnLeftString = leftString(i)
		res@gsnRightString = rightString(i)
		plot(i) = gsn_csm_contour_map_ce(wks,Comp_mean(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.35
	resP@pmLabelBarHeightF = 0.045
	; resP@pmLabelBarParallelPosF = 0.015
	resP@lbLabelFontHeightF = 0.015
	gsn_panel(wks,plot,(/3,1/),resP)  


