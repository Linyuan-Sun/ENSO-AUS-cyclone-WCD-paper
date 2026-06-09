  season = "SON"                                                       ; 定义季节为秋季（SON）
  lonMin = 0                                                         ; 经度范围最小值：-180°
  lonMax = 360                                                          ; 经度范围最大值：180°
  latMin = -90                                                          ; 纬度范围最小值：-90°
  latMax = 0    

f = addfile("climatology_"+season+".nc","r")  ; 读取已计算好的轨迹密度文件
track_density = lonFlip(f->track_density)                     ; 获取轨迹密度变量

  ;==== 绘制 ====
  wks = gsn_open_wks("pdf","Fig6_CTD")               ; 打开绘图工作站，输出 PDF 文件
  res = True                                                            ; 创建绘图资源对象
  res@cnFillOn = True                                                   ; 开启填充颜色（绘制等高线图）
  res@cnLinesOn = False                                                 ; 关闭等高线
;   res@gsnAddCyclic = False                                              ; 禁用循环经度（数据非全球周期性）                                                ; 设置地图经度最大值
  res@cnLevelSelectionMode = "ManualLevels"                             ; 手动设置色标范围
  res@cnMinLevelValF = 0.5       ; 最小值
  res@cnMaxLevelValF = 5      ; 最大值
  res@cnLevelSpacingF = 0.5     ; 等值线间隔
  res@cnFillMode = "RasterFill"
  res@cnFillPalette = "WhiteYellowOrangeRed"                                   ; 设置颜色表
; 设置地图投影和经纬度范围
res@mpMinLatF = latMin              ; 最小纬度（北极区域）
res@mpMaxLatF = latMax              ; 最大纬度（北极）
res@mpMinLonF = lonMin            ; 最小经度
res@mpMaxLonF = lonMax             ; 最大经度
res@mpCenterLonF = 180
res@gsnMajorLonSpacing = 60

res@tmXTOn = False            ; 关闭X轴标签
res@tmYROn = False            ; 关闭Y轴标签
res@pmLabelBarHeightF = 0.05 ; 设置标签栏高度
res@pmLabelBarWidthF = 0.6  ; 设置标签栏宽度
res@pmLabelBarOrthogonalPosF = 0.25
res@pmLabelFontHeightF = 0.012

res@tmXBLabelFontHeightF = 0.015 ; 设置X轴标签字体大小
res@tmYLLabelFontHeightF = 0.015 ; 设置Y轴标签字体大小

res@gsnLeftString = "(b) Cyclone track density" ; 设置子图标签
res@gsnLeftStringFontHeightF = 0.015 ; 设置子图标签字体大小
res@gsnLeftStringOrthogonalPosF = 0.01

res@gsnRightString = "SON climatology 1979-2022" ; 设置子图标签
res@gsnRightStringFontHeightF = 0.015 ; 设置子图标签字体大小
res@gsnRightStringOrthogonalPosF = 0.01

  plot = gsn_csm_contour_map(wks, track_density, res)                   ; 绘制轨迹密度等高线地图