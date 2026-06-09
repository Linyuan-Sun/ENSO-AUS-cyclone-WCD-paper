type = "pattern" ;"pattern" "nonpattern" "season"
pltType = "pdf"
pltName = "Fig11_vq"

siglvl = 0.95

latMin = -60
latMax = 0
lonMin = 100
lonMax = 190

FontHeightF = 0.03

var_EP = "var_EPElNino_"+type+"_mean"
var_CP = "var_CPElNino_"+type+"_mean"
var_La = "var_LaNina_"+type+"_mean"

var_EP_prob = "var_EPElNino_"+type+"_prob"
var_CP_prob = "var_CPElNino_"+type+"_prob"
var_La_prob = "var_LaNina_"+type+"_prob"

;-----div
f = addfile("pattern_out_mvimd.nc","r");
EP_div = f->$var_EP$({latMin:latMax},{lonMin:lonMax})
CP_div = f->$var_CP$({latMin:latMax},{lonMin:lonMax})
La_div = f->$var_La$({latMin:latMax},{lonMin:lonMax})
lat = EP_div&lat
lon = EP_div&lon
nlat = dimsizes(lat)
nlon = dimsizes(lon)
div = new((/3,nlat,nlon/),"float")
div(0,:,:) = EP_div*1e5
div(1,:,:) = CP_div*1e5
div(2,:,:) = La_div*1e5
div!0 = "type"
div!1 = "lat"
div!2 = "lon"
div&lat = lat
div&lon = lon

EP_prob = f->$var_EP_prob$({latMin:latMax},{lonMin:lonMax})
CP_prob = f->$var_CP_prob$({latMin:latMax},{lonMin:lonMax})
La_prob = f->$var_La_prob$({latMin:latMax},{lonMin:lonMax})
Comp_prob = new((/3,nlat,nlon/),"float")
Comp_prob(0,:,:) = EP_prob
Comp_prob(1,:,:) = CP_prob
Comp_prob(2,:,:) = La_prob

;-----uq
f = addfile("pattern_out_viwve.nc","r");
EP_uq = f->$var_EP$({latMin:latMax},{lonMin:lonMax})
CP_uq = f->$var_CP$({latMin:latMax},{lonMin:lonMax})
La_uq = f->$var_La$({latMin:latMax},{lonMin:lonMax})
lat = EP_uq&lat
lon = EP_uq&lon
nlat = dimsizes(lat)
nlon = dimsizes(lon)
uq = new((/3,nlat,nlon/),"float")
uq(0,:,:) = EP_uq
uq(1,:,:) = CP_uq
uq(2,:,:) = La_uq
uq!0 = "type"
uq!1 = "lat"
uq!2 = "lon"
uq&lat = lat
uq&lon = lon

;-----vq
f = addfile("pattern_out_viwvn.nc","r");
EP_vq = f->$var_EP$({latMin:latMax},{lonMin:lonMax})
CP_vq = f->$var_CP$({latMin:latMax},{lonMin:lonMax})
La_vq = f->$var_La$({latMin:latMax},{lonMin:lonMax})
lat = EP_vq&lat
lon = EP_vq&lon
nlat = dimsizes(lat)
nlon = dimsizes(lon)
vq = new((/3,nlat,nlon/),"float")
vq(0,:,:) = EP_vq
vq(1,:,:) = CP_vq
vq(2,:,:) = La_vq
vq!0 = "type"
vq!1 = "lat"
vq!2 = "lon"
vq&lat = lat
vq&lon = lon

threshold = 15
wsd = sqrt(uq*uq + vq*vq)
vq = where(wsd .le. threshold, vq@_FillValue, vq)

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
	res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3/)
	cmap = read_colormap_file("MPL_PuOr")
    res@cnFillPalette = cmap(::-1,:)
	res@tmYLLabelFontHeightF = FontHeightF
	res@tmXBLabelFontHeightF = FontHeightF
	res@gsnMajorLatSpacing = 15
	res@tmXBMode = "Explicit"
	res@tmXBValues = (/100,115,130,145, 160, 175, 190/)
	res@tmXBLabels = (/"100E","115E","130E","145E","160E","175E","170W"/)
	res@tmXBMinorValues = ispan(100,190,5)
	res@mpFillOn               = False
	res@mpOutlineOn           = True
	res@mpOutlineDrawOrder   = "PostDraw"
	
	res@gsnScalarContour     = True
	res@vcRefMagnitudeF           = 60.0             ; define vector ref mag
	res@vcRefLengthF              = 0.06           ; define length of vec ref
	res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
	res@vcRefAnnoFontHeightF	 = 0.02           ; change ref vec font size
	res@vcGlyphStyle            = "CurlyVector"     ; turn on curly vectors
	res@vcLineArrowColor        = "black"           ; change vector color
	res@vcLineArrowThicknessF   = 1.5               ; change vector thickness
	res@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last
	res@vcMinDistanceF = 0.035                 ; min distance between vectors
	res@vcLineArrowHeadMaxSizeF = 0.02
	res@vcLineArrowHeadMinSizeF = 0.005

	res@gsnLeftStringOrthogonalPosF = 0.01
	res@gsnRightStringOrthogonalPosF = 0.01
	res@gsnRightStringFontHeightF = FontHeightF
	res@gsnLeftStringFontHeightF = FontHeightF
	ntilde  = "n~H-13V2F35~D~FV-2H3~"
	rightString = (/"EP El Ni"+ntilde+"o pattern","CP El Ni"+ntilde+"o pattern","La Ni"+ntilde+"a pattern"/)
	leftString = (/"(d) moist. flux","(e) moist. flux","(f) moist. flux"/)
	
	plot = new(3, graphic)
	plot_prob = new(3, graphic)
	do i = 0, 2	
		res@gsnLeftString = leftString(i)
		res@gsnRightString = rightString(i)
		plot(i) = gsn_csm_vector_scalar_map_ce(wks, uq(i,:,:), vq(i,:,:), div(i,:,:), res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.35
	resP@pmLabelBarHeightF = 0.045
	resP@lbLabelFontHeightF = 0.015
	gsn_panel(wks,plot,(/3,1/),resP)  


