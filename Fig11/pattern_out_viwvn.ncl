; 基本思路
; 1. 读取事件索引文件，获取EP El Nino、CP El Nino和La Nina的年份和起始-结束时刻位置
; 2. 生成各事件年的8月17日-12月15日时间序列
; 3. 遍历每个事件年份，读取对应的天气数据文件
; 4. 把所有事件和相应的气旋日数分别合并到一起，计算频率
; 5. 用气旋数据涵盖年份的季节窗口（SON）来计算频率气候态
; 5. 输出结果到新的NetCDF文件

season = "SON"  ; 季节，SON表示9-11月
startday = 901
endday = 1130
season_window = 91
type = "cyclone"  ; cyclone表示气旋，anticyclone表示反气旋

var = "viwvn"

; ENSO events
EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016, 2024/)

CP_ElNino_year = (/1987, 1992, 1995, 2003, 2005, 2010, 2019/)

EP_LaNina_year = (/1996, 2000, 2008, 2021/)
CP_LaNina_year = (/1984, 1989, 1999, 2011, 2012, 2022, 2023/)
LaNina_year = array_append_record(EP_LaNina_year,CP_LaNina_year,0)

; load event index
f = addfile("/g/data/m35/ls3248/ENSO_cyclone/paper/Fig2/event_identify.nc.nc", "r")
nEPElNino = f->nEPElNino
nCPElNino = f->nCPElNino
nLaNina = f->nLaNina
start_date_EPElNino = f->start_date_EPElNino
end_date_EPElNino = f->end_date_EPElNino
start_date_CPElNino = f->start_date_CPElNino
end_date_CPElNino = f->end_date_CPElNino
start_date_LaNina = f->start_date_LaNina
end_date_LaNina = f->end_date_LaNina
start_idx_EPElNino = f->start_idx_EPElNino
end_idx_EPElNino = f->end_idx_EPElNino
start_idx_CPElNino = f->start_idx_CPElNino
end_idx_CPElNino = f->end_idx_CPElNino
start_idx_LaNina = f->start_idx_LaNina
end_idx_LaNina = f->end_idx_LaNina

;-----计算各事件和非事件窗口的总天数
EPElNino_event_total_days = sum(end_idx_EPElNino - start_idx_EPElNino + 1) 
CPElNino_event_total_days = sum(end_idx_CPElNino - start_idx_CPElNino + 1)
LaNina_event_total_days = sum(end_idx_LaNina - start_idx_LaNina + 1)

EPElNino_nonevent_total_days = season_window*dimsizes(EP_ElNino_year) - EPElNino_event_total_days
CPElNino_nonevent_total_days = season_window*dimsizes(CP_ElNino_year) - CPElNino_event_total_days
LaNina_nonevent_total_days = season_window*dimsizes(LaNina_year) - LaNina_event_total_days

; print("EPElNino_event_total_days: "+EPElNino_event_total_days)
; print("CPElNino_event_total_days: "+CPElNino_event_total_days)
; print("LaNina_event_total_days: "+LaNina_event_total_days)
; print("EPElNino_nonevent_total_days: "+EPElNino_nonevent_total_days)
; print("CPElNino_nonevent_total_days: "+CPElNino_nonevent_total_days)
; print("LaNina_nonevent_total_days: "+LaNina_nonevent_total_days)
; exit

;引入变量数据
nlat = 180
nlon = 360
year_start = 1979
year_end = 2024
total_year = ispan(year_start, year_end, 1)
ntotal_year = dimsizes(total_year)

;-----EPElNino
print("EPElNino:-----")
;型事件期间内的变量和
var_EPElNino = new((/nEPElNino,nlat,nlon/),"float")
do i = 0,nEPElNino-1
    year = start_date_EPElNino(i)/10000
    f = addfile("/g/data/if69/ls3248/data/ERA5/daily/"+var+"/raw/"+var+"."+season+".anom.raw."+year+".nc","r")
    time = cd_calendar(f->time, 2)
    ind_start = ind(time.eq.start_date_EPElNino(i))
    ind_end = ind(time.eq.end_date_EPElNino(i))
    var_EPElNino(i,:,:) = dim_sum_n_Wrap(f->$var$(ind_start:ind_end,:,:),0) ;同一年内各事件相加
end do
var_EPElNino_total = dim_sum_n_Wrap(var_EPElNino,0) ;不同年份再相加
var_EPElNino_pattern_mean = var_EPElNino_total/EPElNino_event_total_days ;型事件变量均值
copy_VarCoords(var_EPElNino_total, var_EPElNino_pattern_mean)
delete(time)

;-----CPElNino
print("CPElNino:-----")
;型事件期间内的变量和
var_CPElNino = new((/nCPElNino,nlat,nlon/),"float")
do i = 0,nCPElNino-1
    year = start_date_CPElNino(i)/10000
    f = addfile("/g/data/if69/ls3248/data/ERA5/daily/"+var+"/raw/"+var+"."+season+".anom.raw."+year+".nc","r")
    time = cd_calendar(f->time, 2)
    ind_start = ind(time.eq.start_date_CPElNino(i))
    ind_end = ind(time.eq.end_date_CPElNino(i))
    var_CPElNino(i,:,:) = dim_sum_n_Wrap(f->$var$(ind_start:ind_end,:,:),0) ;同一年内各事件相加
end do
var_CPElNino_total = dim_sum_n_Wrap(var_CPElNino,0) ;不同年份再相加
var_CPElNino_pattern_mean = var_CPElNino_total/CPElNino_event_total_days ;型事件变量均值
copy_VarCoords(var_CPElNino_total, var_CPElNino_pattern_mean)
delete(time)

;-----LaNina
print("LaNina:-----")
;型事件期间内的变量和
var_LaNina = new((/nLaNina,nlat,nlon/),"float")
do i = 0,nLaNina-1
    year = start_date_LaNina(i)/10000
    f = addfile("/g/data/if69/ls3248/data/ERA5/daily/"+var+"/raw/"+var+"."+season+".anom.raw."+year+".nc","r")
    time = cd_calendar(f->time, 2)
    ind_start = ind(time.eq.start_date_LaNina(i))
    ind_end = ind(time.eq.end_date_LaNina(i))
    var_LaNina(i,:,:) = dim_sum_n_Wrap(f->$var$(ind_start:ind_end,:,:),0) ;同一年内各事件相加
end do
var_LaNina_total = dim_sum_n_Wrap(var_LaNina,0) ;不同年份再相加
var_LaNina_pattern_mean = var_LaNina_total/LaNina_event_total_days ;型事件变量均值
copy_VarCoords(var_LaNina_total, var_LaNina_pattern_mean)
delete(time)

;-----output the Comp_avg_hgt
out_file = "pattern_out_viwvn.nc"

setfileoption("nc","Format","LargeFile")

if (fileexists(out_file)) then
    system("rm -f "+out_file)
end if
f_out = addfile(out_file,"c")
f_out->var_EPElNino_pattern_mean = var_EPElNino_pattern_mean
f_out->var_CPElNino_pattern_mean = var_CPElNino_pattern_mean
f_out->var_LaNina_pattern_mean = var_LaNina_pattern_mean