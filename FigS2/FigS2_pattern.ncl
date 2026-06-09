season = "SON"                                                       ; 定义季节为秋季（SON）
lonMin = 0                                                         ; 经度范围最小值：-180°
lonMax = 360                                                          ; 经度范围最大值：180°
latMin = -90                                                          ; 纬度范围最小值：-90°
latMax = 0    

; 气候态
f = addfile("/g/data/m35/ls3248/ENSO_cyclone/synoptic_system/track_density/climatology/climatology_SON.nc","r")  ; 读取已计算好的轨迹密度文件
track_density_clim = lonFlip(f->track_density)                     ; 获取轨迹密度变量

; ENSO
f = addfile("/g/data/m35/ls3248/ENSO_cyclone/synoptic_system/track_density/event/track_density_pattern.nc","r")  ; 读取已计算好的轨迹密度文件
track_density_EPElNino = lonFlip(f->track_density_EPElNino)
track_density_CPElNino = lonFlip(f->track_density_CPElNino)
track_density_LaNina = lonFlip(f->track_density_LaNina)

; 计算异常值
anom_EPElNino = track_density_EPElNino - track_density_clim
anom_CPElNino = track_density_CPElNino - track_density_clim
anom_LaNina = track_density_LaNina - track_density_clim
copy_VarCoords(track_density_clim, anom_EPElNino)
copy_VarCoords(track_density_clim, anom_CPElNino)
copy_VarCoords(track_density_clim, anom_LaNina)

;merge
lat = track_density_clim&lat
lon = track_density_clim&lon
nlat = dimsizes(lat)
nlon = dimsizes(lon)
anom = new((/3,nlat,nlon/),"float")
anom(0,:,:) = anom_EPElNino
anom(1,:,:) = anom_CPElNino
anom(2,:,:) = anom_LaNina
anom!1 = "lat"
anom!2 = "lon"
anom&lat = lat
anom&lon = lon

; 绘图
wks = gsn_open_wks("pdf","FigS2_pattern")

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
res@cnFillMode = "RasterFill"

res@mpMinLatF = latMin              ; 最小纬度（北极区域）
res@mpMaxLatF = latMax              ; 最大纬度（北极）
res@mpMinLonF = lonMin            ; 最小经度
res@mpMaxLonF = lonMax             ; 最大经度
res@mpCenterLonF = (lonMin+lonMax)/2
res@mpLandFillColor = "white" ; 设置地图填充颜色
; res@mpFillDrawOrder = "PostDraw" ; 设置填充绘制顺序
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/-6,-5,-4,-3,-2,-1,1,2,3,4,5,6/)/2.
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
LeftString = (/"(d) Cyclone track density", "(e) Cyclone track density", "(f) Cyclone track density"/)
ntilde  = "n~H-13V2F35~D~FV-2H3~"
RightString = (/"EP El Ni"+ntilde+"o pattern", "CP El Ni"+ntilde+"o pattern", "La Ni"+ntilde+"a pattern"/)
do i = 0, 2
		res@gsnLeftString = LeftString(i)
		res@gsnRightString = RightString(i)
  plot(i) = gsn_csm_contour_map_ce(wks,anom(i,:,:),res)
end do
resP                  = True   
resP@gsnMaximize = True
resP@gsnPanelLabelBar = True   
resP@pmLabelBarWidthF = 1.
resP@pmLabelBarHeightF = 0.07
resP@lbLabelFontHeightF = 0.02
gsn_panel(wks,plot,(/3,1/),resP)  