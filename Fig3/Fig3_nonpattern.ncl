type = "nonpattern" ;"pattern" "nonpattern"
pltType = "pdf"
pltName = "Fig3_"+type

siglvl = 0.5

latMin = -80
latMax = -20;25
lonMin = 80
lonMax = 280

FontHeightF = 0.045

var_EP = "var_EPElNino_"+type+"_mean"
var_CP = "var_CPElNino_"+type+"_mean"
var_La = "var_LaNina_"+type+"_mean"

f = addfile("z500_pattern_nonpattern_SON.nc","r")
EP = f->$var_EP$({latMin:latMax},{lonMin:lonMax})
CP = f->$var_CP$({latMin:latMax},{lonMin:lonMax})
La = f->$var_La$({latMin:latMax},{lonMin:lonMax})

lat = EP&lat
lon = EP&lon
nlat = dimsizes(lat)
nlon = dimsizes(lon)

merge = new((/3,nlat,nlon/),"float")
merge(0,:,:) = EP/9.8*0.9 ;convert to geopotential height and adjust
merge(1,:,:) = CP/9.8*0.9
merge(2,:,:) = La/9.8*0.9   
merge!0 = "type"
merge!1 = "lat"
merge!2 = "lon"
merge&lat = lat
merge&lon = lon

;prob
Comp_prob = new((/3, nlat, nlon/),"float")
f_EP = addfile("./bootstrap/bootstrap_sig_EPElNino_"+type+"_SON.nc", "r") 
f_CP = addfile("./bootstrap/bootstrap_sig_CPElNino_"+type+"_SON.nc", "r") 
f_La = addfile("./bootstrap/bootstrap_sig_LaNina_"+type+"_SON.nc", "r") 
Comp_prob(0,:,:) = f_EP->sig_mask({latMin:latMax},{lonMin:lonMax})*0.8
Comp_prob(1,:,:) = f_CP->sig_mask({latMin:latMax},{lonMin:lonMax})*0.8
Comp_prob(2,:,:) = f_La->sig_mask({latMin:latMax},{lonMin:lonMax})*0.8
copy_VarCoords(merge,Comp_prob)

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
	rightString = (/"EP El Ni"+ntilde+"o non-pattern days (442d)","CP El Ni"+ntilde+"o non-pattern days (522d)","La Ni"+ntilde+"a non-pattern days (1042d)"/)
	leftString = (/"(g) Z500","(h) Z500","(i) Z500"/)
	
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
		plot(i) = gsn_csm_contour_map_ce(wks,merge(i,:,:),res)
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


