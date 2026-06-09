FontHeightF = 0.04

f = addfile("decomp_ENSO.nc","r")
pc_EPElNino = f->pc_EPElNino
pnc_EPElNino = f->pnc_EPElNino
cyclone_days_num_EPElNino = f->cyclone_days_num_EPElNino
noncyclone_days_num_EPElNino = f->noncyclone_days_num_EPElNino

pc_CPElNino = f->pc_CPElNino
pnc_CPElNino = f->pnc_CPElNino
cyclone_days_num_CPElNino = f->cyclone_days_num_CPElNino
noncyclone_days_num_CPElNino = f->noncyclone_days_num_CPElNino	

pc_LaNina = f->pc_LaNina
pnc_LaNina = f->pnc_LaNina
cyclone_days_num_LaNina = f->cyclone_days_num_LaNina
noncyclone_days_num_LaNina = f->noncyclone_days_num_LaNina


f = addfile("decomp_clim.nc","r")
pc_clim = f->pc_clim
pnc_clim = f->pnc_clim

cyclone_EPElNino = pc_EPElNino - pc_clim
noncyclone_EPElNino = pnc_EPElNino - pnc_clim
cyclone_CPElNino = pc_CPElNino - pc_clim
noncyclone_CPElNino = pnc_CPElNino - pnc_clim
cyclone_LaNina = pc_LaNina - pc_clim
noncyclone_LaNina = pnc_LaNina - pnc_clim

copy_VarCoords(pc_EPElNino, cyclone_EPElNino)
copy_VarCoords(pc_EPElNino, noncyclone_EPElNino)
copy_VarCoords(pc_EPElNino, cyclone_CPElNino)
copy_VarCoords(pc_EPElNino, noncyclone_CPElNino)
copy_VarCoords(pc_EPElNino, cyclone_LaNina)
copy_VarCoords(pc_EPElNino, noncyclone_LaNina)

EP = cyclone_EPElNino + noncyclone_EPElNino
CP = cyclone_CPElNino + noncyclone_CPElNino
La = cyclone_LaNina + noncyclone_LaNina
copy_VarCoords(pc_EPElNino, EP)
copy_VarCoords(pc_EPElNino, CP)
copy_VarCoords(pc_EPElNino, La)

; ;pattern event precip
; f = addfile("/g/data/m35/ls3248/ENSO_synop/comp_pattern_based/pattern_check/precip/precip_pattern_nonpattern_SON.nc","r")
; EP = f->var_EPElNino_pattern_mean
; CP = f->var_CPElNino_pattern_mean
; La = f->var_LaNina_pattern_mean

nlat = dimsizes(pc_EPElNino(:,0))
nlon = dimsizes(pc_EPElNino(0,:))
term = new((/9,nlat,nlon/),"float")
term(0,:,:) = EP
term(1,:,:) = cyclone_EPElNino*1.15
term(2,:,:) = noncyclone_EPElNino*0.85
term(3,:,:) = CP
term(4,:,:) = cyclone_CPElNino*1.15
term(5,:,:) = noncyclone_CPElNino*0.85
term(6,:,:) = La
term(7,:,:) = cyclone_LaNina*1.15
term(8,:,:) = noncyclone_LaNina*0.85

term!1 = "lat"
term!2 = "lon"
term&lat = pc_EPElNino&lat
term&lon = pc_EPElNino&lon

latMin = -45
latMax = -10
lonMin = 135
lonMax = 155

wks = gsn_open_wks("pdf","relative_contribution")  ; open a workstation
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

	res@mpMinLatF = latMin 
    res@mpMaxLatF = latMax
    res@mpMinLonF = lonMin   
    res@mpMaxLonF = lonMax
    res@mpCenterLonF = (lonMin+lonMax)/2.
	res@mpShapeMode = "FreeAspect"
	res@vpHeightF         = 0.8                    ; change aspect ratio of plot
	res@vpWidthF          = 0.6
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
	res@gsnMajorLonSpacing = 10

	res@gsnLeftStringOrthogonalPosF = 0.01
	res@gsnRightStringOrthogonalPosF = 0.01
	res@gsnRightStringFontHeightF = FontHeightF
	res@gsnLeftStringFontHeightF = FontHeightF
	ntilde  = "n~H-13V2F35~D~FV-2H3~"
	rightString = (/"EP El Ni"+ntilde+"o","EP El Ni"+ntilde+"o", "EP El Ni"+ntilde+"o", \
					"CP El Ni"+ntilde+"o","CP El Ni"+ntilde+"o","CP El Ni"+ntilde+"o", \
					"La Ni"+ntilde+"a","La Ni"+ntilde+"a","La Ni"+ntilde+"a"/)
	leftString = (/"(a) P~B~anom~","(a) P~S1~c~B4~anom~","(d) P~S2~nc~B4~anom~", \
					"(b) P~B~anom~","(b) P~S1~c~B4~anom~","(e) P~S2~nc~B4~anom~", \
					"(c) P~B~anom~","(c) P~S1~c~B4~anom~","(f) P~S2~nc~B4~anom~"/)

	YLLabel = (/True,False,False,True,False,False,True,False,False/)
	XBLabel = (/False,False,False,False,False,False,True,True,True/)
	plot = new(9, graphic)
    do i = 0, 8
		res@tmYLLabelsOn = YLLabel(i)
		res@tmXBLabelsOn = XBLabel(i)
        res@gsnLeftString = leftString(i)
		res@gsnRightString = rightString(i)
		plot(i) = gsn_csm_contour_map_ce(wks,term(i,:,:),res)
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.6
	resP@pmLabelBarHeightF = 0.045
	; resP@pmLabelBarParallelPosF = 0.015
	resP@pmLabelBarOrthogonalPosF = -0.03
	resP@lbLabelFontHeightF = 0.015
	gsn_panel(wks,plot,(/3,3/),resP)  


