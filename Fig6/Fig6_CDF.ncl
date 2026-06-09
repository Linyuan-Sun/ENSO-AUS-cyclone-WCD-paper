sphere = "SH"
latMax = 0
latMin = -90;-60
lonMax = 360;200
lonMin = 0;90

season = "SON"
mon1 = 9
mon2 = 10
mon3 = 11
days_window = 91
threshold = 6

cyclone_data_year = ispan(1970, 2022, 1)

year_start = 1979
year_end = 2022
total_year = ispan(year_start, year_end, 1)
type = "cyclone"
pltTpye = "pdf"
pltName = "Fig6_CDF"

;每月气旋日总数
f = addfile("/g/data/if69/ls3248/data/weather/"+type+"/"+season+"/mon/"+type+"."+season+".mon.nc", "r") ;total days in each month
time = cd_calendar(f->time, 0)
;提取月份
ind_month = ind(time(:,1).eq.mon1 .or. time(:,1).eq.mon2 .or. time(:,1).eq.mon3) 
flag_temp_temp = f->flag(ind_month, :, :)
delete([/time,ind_month/])


;转换为年-月-经纬度形状
flag_reshape = reshape(flag_temp_temp, (/dimsizes(flag_temp_temp(:,0,0))/3, 3, dimsizes(flag_temp_temp(0,:,0)), dimsizes(flag_temp_temp(0,0,:))/))
flag_reshape!0 = "year"
flag_reshape&year = cyclone_data_year

flag = flag_reshape({year_start:year_end},:,:,:)*1. ;提取实际气候态年份，转换为浮点型

flag_season = dim_sum_n(flag, 1)/days_window*100. ; 计算季节内的气旋日频率（百分比）
flag_season!0 = "time" ; 设置时间维度名称
flag_season!1 = "lat" ; 设置纬度维度名称
flag_season!2 = "lon" ; 设置经度维度名称
flag_season&time = total_year
flag_season&lat = flag_temp_temp&lat
flag_season&lon = flag_temp_temp&lon
flag_season_clim = dim_avg_n_Wrap(flag_season, 0)
; flag_season_clim = where(flag_season_clim.lt.1, flag_season_clim@_FillValue, flag_season_clim) ; 频率不能小于阈值
; printVarSummary(flag_season_clim)

wks = gsn_open_wks(pltTpye, pltName)
; 设置绘图资源
res = True
res@gsnFrame = False            ; 不立即绘制
res@gsnDraw  = False            ; 不立即显示
res@cnFillOn = True             ; 开启填色
res@cnLinesOn = False           ; 关闭等高线（仅显示填色）
res@cnFillPalette = "WhiteYellowOrangeRed";"MPL_YlOrRd" ; 选择颜色表
res@cnLevelSelectionMode = "ExplicitLevels" ; 手动设置等值线级别
res@cnLevels = (/1, 3,5,7,9,11,13,15,20, 30,50/) ; 设置等值线级别

; 设置地图投影和经纬度范围
res@mpMinLatF = latMin              ; 最小纬度（北极区域）
res@mpMaxLatF = latMax              ; 最大纬度（北极）
res@mpMinLonF = lonMin            ; 最小经度
res@mpMaxLonF = lonMax             ; 最大经度
res@mpCenterLonF = (lonMin+lonMax)/2
res@gsnMajorLonSpacing = 60

res@tmXTOn = False            ; 关闭X轴标签
res@tmYROn = False            ; 关闭Y轴标签
res@pmLabelBarHeightF = 0.05 ; 设置标签栏高度
res@pmLabelBarWidthF = 0.6  ; 设置标签栏宽度
res@pmLabelBarOrthogonalPosF = 0.25
res@pmLabelFontHeightF = 0.012

res@tmXBLabelFontHeightF = 0.015 ; 设置X轴标签字体大小
res@tmYLLabelFontHeightF = 0.015 ; 设置Y轴标签字体大小

res@gsnLeftString = "(a) Cyclone day frequency (CDF)" ; 设置子图标签
res@gsnLeftStringFontHeightF = 0.015 ; 设置子图标签字体大小
res@gsnLeftStringOrthogonalPosF = 0.01

res@gsnRightString = "SON climatology 1979-2022" ; 设置子图标签
res@gsnRightStringFontHeightF = 0.015 ; 设置子图标签字体大小
res@gsnRightStringOrthogonalPosF = 0.01


; 绘制填色图
plot = gsn_csm_contour_map(wks, flag_season_clim, res)

; 显示图形
draw(plot)
frame(wks)