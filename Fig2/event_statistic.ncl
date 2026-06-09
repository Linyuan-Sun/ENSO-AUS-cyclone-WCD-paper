;-----SON define
; EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016, 2024/)
; CP_ElNino_year = (/1987, 1992, 1995, 2003, 2005, 2010, 2019/)
; ; EP_LaNina_year = (/1985, 1996, 2000, 2008, 2018, 2021/)
; LaNina_year = (/1984, 1989, 1999, 2001, 2011, 2012, 2022, 2023/)

EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016/)
CP_ElNino_year = (/1987, 1992, 1995, 2003, 2010, 2019/)
EP_LaNina_year = (/1996, 2000, 2008, 2021/)
CP_LaNina_year = (/1984, 1989, 1999, 2011, 2012, 2022, 2023/)
LaNina_year = array_append_record(EP_LaNina_year,CP_LaNina_year,0)
ElNino_year = array_append_record(EP_ElNino_year, CP_ElNino_year, 0)
ENSO_year = array_append_record(ElNino_year, LaNina_year, 0)
qsort(ENSO_year)

f = addfile("event_identify.nc", "r")
event_year = f->year
max_idx = f->max_idx
start_idx = f->start_idx
end_idx = f->end_idx

; 构造y轴
y_axis = ispan(1, dimsizes(ENSO_year), 1)

; 构造起始点、中间最值点数组
event_year_num = dimsizes(event_year)
x = new((/event_year_num,3/),"float")
y = x
do i = 0,event_year_num-1
    x(i,0) = (/start_idx(i)/)+1
    x(i,1) = (/max_idx(i)/)+1
    x(i,2) = (/end_idx(i)/)+1
    y(i,:) = y_axis(ind(ENSO_year-1.eq.event_year(i)))
end do

; 创建工作站
wks_type = "pdf"
wks = gsn_open_wks(wks_type, "event_statistic")

; 动态设置y轴标签
tmYLValues = y_axis
tmYLLabels = y_axis+""

ntilde  = "n~H-13V2F35~D~FV-2H3~"

do i = 0,dimsizes(ENSO_year)-1
    if any(EP_ElNino_year.eq.ENSO_year(i)) then
        tmYLLabels(i) = ENSO_year(i)-1+" EP El Ni"+ntilde+"o"
    elseif any(CP_ElNino_year.eq.ENSO_year(i))
        tmYLLabels(i) = ENSO_year(i)-1+" CP El Ni"+ntilde+"o"
    elseif any(LaNina_year.eq.ENSO_year(i))
        tmYLLabels(i) = ENSO_year(i)-1+" La Ni"+ntilde+"a"
    end if
end do

; 设置绘图资源
res = True
res@gsnFrame = False
res@gsnDraw = False
res@gsnMaximize = True
res@trYMinF = 0       ; y轴最小值
res@trYMaxF = dimsizes(ENSO_year)+1       ; y轴最大值
res@trXMinF = 0          ; x轴最小值
res@trXMaxF = 120        ; x轴最大值
res@tmYLMode = "Explicit" ; 显式设置y轴标签
res@tmYLValues = tmYLValues ; y轴刻度值
res@tmYLLabels = tmYLLabels ; y轴标签
res@tiYAxisString = "ENSO events"
res@tiYAxisFontHeightF = 0.025
res@tmXTOn = False
res@tmYROn = False
res@vpWidthF = 1.2
res@vpHeightF = 1;0.8

; 添加网格线
res@tmXMajorGrid = True              ; 启用x轴主要网格线
res@tmYMajorGrid = True             ; 启用y轴主要网格线
res@tmXMajorGridLineColor = "gray"   ; x轴网格线颜色
res@tmYMajorGridLineColor = "gray"  ; y轴网格线颜色
res@tmXMajorGridLineDashPattern = 2  ; x轴网格线为虚线
res@tmYMajorGridLineDashPattern = 2 ; y轴网格线为虚线
res@tmXMajorGridThicknessF = 2     ; x轴网格线粗细
res@tmYMajorGridThicknessF = 2    ; y轴网格线粗细
res@tmGridDrawOrder = "PreDraw"      ; 网格线绘制在数据之前

; Set x-axis (days: 1 to 121)
res@trXMinF = 0
res@trXMaxF = 122

major_ticks = (/16, 46, 77, 107/) 
major_labels = (/"1st Sep", "1st Oct", "1st Nov", "1st Dec"/)
res@tmXBMode = "Explicit"                  ; Set explicit tick marks
res@tmXBValues = major_ticks               ; Major tick locations
res@tmXBLabels = major_labels              ; Major tick labels
res@tmXBMinorOn = True                     ; Enable minor ticks
res@tmXBMinorValues = ispan(1,121,1) 

res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

; 创建空白画布
plot = gsn_csm_xy(wks, (/0, 0/), (/0, 0/), res)

; 线段资源
line_res = True
line_res@gsLineThicknessF = 5.0

; 关键点资源
marker_res = True
marker_res@gsMarkerIndex = 16 ; 圆形标记
marker_res@gsMarkerSizeF = 0.015

draw(plot) ;先 draw plot（网格线在此时画出）

; 再循环绘制线段和关键点（覆盖在网格线之上）
do i = 0,event_year_num-1
    if (any(ENSO_year.eq.event_year(i)+1)) then
        ; 根据年份设置颜色
        if (any(EP_ElNino_year.eq.event_year(i)+1)) then
            line_res@gsLineColor = "red"
            marker_res@gsMarkerColor = "red"
        elseif any(CP_ElNino_year.eq.event_year(i)+1)
            line_res@gsLineColor = "orange"
            marker_res@gsMarkerColor = "orange"
        elseif any(LaNina_year.eq.event_year(i)+1)
            line_res@gsLineColor = "blue"
            marker_res@gsMarkerColor = "blue"
        end if
        ; 绘制线段
        gsn_polyline(wks, plot, x(i,:), y(i,:), line_res)
        ; 绘制关键点（中间点）
        gsn_polymarker(wks, plot, x(i,1), y(i,1), marker_res)
    end if
end do

frame(wks)