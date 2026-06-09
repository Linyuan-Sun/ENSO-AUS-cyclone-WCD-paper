; 基本思路
; 1. 读取事件索引文件，获取EP El Nino、CP El Nino和La Nina各事件的起止日期及索引位置
; 2. 计算各ENSO类型的pattern事件总天数
; 3. 遍历各ENSO类型的每个事件年份，读取对应ERA5逐日变量文件：
;    a. 提取事件期间单层数据，同步存入逐日数组（用于显著性检验）并累加年际求和
; 4. 计算pattern_mean：事件期间所有日数据的总和除以pattern总天数
; 5. 对pattern_mean进行考虑滞后1自相关的有效自由度单样本t检验
; 6. 将pattern_mean及对应显著性概率场输出到新的NetCDF文件

season = "SON"  ; 季节，SON表示9-11月
startday = 901
endday = 1130
season_window = 91
type = "cyclone"  ; cyclone表示气旋，anticyclone表示反气旋

var = "w"
level = 500

; 空间范围限制（仅南半球，节省内存）
lat_s = -90.
lat_n =   0.
lon_w =   0.
lon_e = 360.

; ENSO events
EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016, 2024/)
CP_ElNino_year = (/1987, 1992, 1995, 2003, 2005, 2010, 2019/)
EP_LaNina_year = (/1996, 2000, 2008, 2021/)
CP_LaNina_year = (/1984, 1989, 1999, 2011, 2012, 2022, 2023/)
LaNina_year = array_append_record(EP_LaNina_year,CP_LaNina_year,0)

; load event index
f = addfile("/g/data/m35/ls3248/ENSO_cyclone/paper/Fig2/event_identify.nc", "r")
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

;-----计算各事件总天数
EPElNino_event_total_days = sum(end_idx_EPElNino - start_idx_EPElNino + 1)
CPElNino_event_total_days = sum(end_idx_CPElNino - start_idx_CPElNino + 1)
LaNina_event_total_days   = sum(end_idx_LaNina   - start_idx_LaNina   + 1)

; 根据空间限制确定实际格点数（从第一个文件中读取）
f_tmp = addfile("/g/data/if69/ls3248/data/ERA5/daily/"+var+"/raw/"+var+"."+season+".anom.raw.1983.nc","r")
lat_tmp = f_tmp->lat({lat_s:lat_n})
lon_tmp = f_tmp->lon({lon_w:lon_e})
nlat = dimsizes(lat_tmp)
nlon = dimsizes(lon_tmp)
delete([/f_tmp, lat_tmp, lon_tmp/])

year_start = 1979
year_end   = 2024
total_year = ispan(year_start, year_end, 1)
ntotal_year = dimsizes(total_year)

; ===========================================================
; 辅助函数：利用滞后1自相关修正有效自由度
; 有效自由度 neff = n * (1 - r1) / (1 + r1)
; 利用esacr_n向量化运算，返回 (nlat,nlon) 的有效自由度场
; ===========================================================
undef("calc_eff_df")
function calc_eff_df(data[*][*][*]:float)  ; data维度: (nday, nlat, nlon)
local nday, r1, neff
begin
    nday = dimsizes(data(:,0,0))
    ; esacr_n 沿第0维（时间）计算滞后自相关，返回 (mxlag+1, nlat, nlon)
    ; 取 [1,:,:] 即滞后1的自相关系数场
    r1_temp = esacr_n(data, 1, 0)   ; r1(2, nlat, nlon)
    r1 = r1_temp(1,:,:)              ; r1(nlat, nlon)
    ; 限制 r1 范围，防止除零或 neff 为负
    r1 = where(r1 .ge.  1.0,  0.999, r1)
    r1 = where(r1 .le. -1.0, -0.999, r1)
    ; 有效自由度：neff = n * (1 - r1) / (1 + r1)
    neff = nday * (1.0 - r1) / (1.0 + r1)
    ; 保证 neff 至少为 2（t 检验最低要求）
    neff = where(neff .lt. 2.0, 2.0, neff)
    return(neff)
end

; ===========================================================
; 辅助函数：单样本 t 检验显著性（双尾）
; 输入: data(nday,nlat,nlon), 返回 prob(nlat,nlon)
; ===========================================================
undef("ttest_1sample_lag1")
function ttest_1sample_lag1(data[*][*][*]:float)
local avg, std, neff, tval, df, t2, a, b, prob
begin
    avg  = dim_avg_n_Wrap(data, 0)
    std  = dim_stddev_n_Wrap(data, 0)
    neff = calc_eff_df(data)
    tval = avg / (std / sqrt(neff))
    df   = neff - 1.0
    df   = where(df .lt. 1.0, 1.0, df)
    t2   = tval^2
    ; 构造与 df 同维度的 a、b 数组，满足 betainc 的维度要求
    a    = df / 2.0
    b    = new(dimsizes(df), float)
    b    = 0.5
    prob = 1.0 - betainc(df/(df+t2), a, b)
    copy_VarCoords(avg, prob)
    return(prob)
end

; ===========================================================
; -----EP El Nino
; ===========================================================
print("EPElNino:-----")

nday_pat_EP = EPElNino_event_total_days
var_EPElNino_patdays = new((/nday_pat_EP,nlat,nlon/),"float")
idx_pat = 0
var_EPElNino = new((/nEPElNino,nlat,nlon/),"float")
do i = 0,nEPElNino-1
    year = start_date_EPElNino(i)/10000
    f = addfile("/g/data/if69/ls3248/data/ERA5/daily/"+var+"/raw/"+var+"."+season+".anom.raw."+year+".nc","r")
    time = cd_calendar(f->time, 2)
    ind_start = ind(time.eq.start_date_EPElNino(i))
    ind_end   = ind(time.eq.end_date_EPElNino(i))
    ndays_i   = ind_end - ind_start + 1
    ; 单层数据直接得到 (nday,nlat,nlon)，无需垂直平均
    var_EPElNino_patdays(idx_pat:idx_pat+ndays_i-1,:,:) = \
        f->$var$(ind_start:ind_end,{level},{lat_s:lat_n},{lon_w:lon_e})
    var_EPElNino(i,:,:) = \
        dim_sum_n_Wrap(f->$var$(ind_start:ind_end,{level},{lat_s:lat_n},{lon_w:lon_e}),0)
    idx_pat = idx_pat + ndays_i
    delete([/time, ind_start, ind_end, ndays_i/])
end do
var_EPElNino_total        = dim_sum_n_Wrap(var_EPElNino,0)
var_EPElNino_pattern_mean = var_EPElNino_total/EPElNino_event_total_days
copy_VarCoords(var_EPElNino_total, var_EPElNino_pattern_mean)
var_EPElNino_pattern_prob = ttest_1sample_lag1(var_EPElNino_patdays)
; 清理EP El Nino中间变量
delete([/var_EPElNino, var_EPElNino_patdays, var_EPElNino_total/])

; ===========================================================
; -----CP El Nino
; ===========================================================
print("CPElNino:-----")

nday_pat_CP = CPElNino_event_total_days
var_CPElNino_patdays = new((/nday_pat_CP,nlat,nlon/),"float")
idx_pat = 0
var_CPElNino = new((/nCPElNino,nlat,nlon/),"float")
do i = 0,nCPElNino-1
    year = start_date_CPElNino(i)/10000
    f = addfile("/g/data/if69/ls3248/data/ERA5/daily/"+var+"/raw/"+var+"."+season+".anom.raw."+year+".nc","r")
    time = cd_calendar(f->time, 2)
    ind_start = ind(time.eq.start_date_CPElNino(i))
    ind_end   = ind(time.eq.end_date_CPElNino(i))
    ndays_i   = ind_end - ind_start + 1
    var_CPElNino_patdays(idx_pat:idx_pat+ndays_i-1,:,:) = \
        f->$var$(ind_start:ind_end,{level},{lat_s:lat_n},{lon_w:lon_e})
    var_CPElNino(i,:,:) = \
        dim_sum_n_Wrap(f->$var$(ind_start:ind_end,{level},{lat_s:lat_n},{lon_w:lon_e}),0)
    idx_pat = idx_pat + ndays_i
    delete([/time, ind_start, ind_end, ndays_i/])
end do
var_CPElNino_total        = dim_sum_n_Wrap(var_CPElNino,0)
var_CPElNino_pattern_mean = var_CPElNino_total/CPElNino_event_total_days
copy_VarCoords(var_CPElNino_total, var_CPElNino_pattern_mean)
var_CPElNino_pattern_prob = ttest_1sample_lag1(var_CPElNino_patdays)
; 清理CP El Nino中间变量
delete([/var_CPElNino, var_CPElNino_patdays, var_CPElNino_total/])

; ===========================================================
; -----La Nina
; ===========================================================
print("LaNina:-----")

nday_pat_LN = LaNina_event_total_days
var_LaNina_patdays = new((/nday_pat_LN,nlat,nlon/),"float")
idx_pat = 0
var_LaNina = new((/nLaNina,nlat,nlon/),"float")
do i = 0,nLaNina-1
    year = start_date_LaNina(i)/10000
    f = addfile("/g/data/if69/ls3248/data/ERA5/daily/"+var+"/raw/"+var+"."+season+".anom.raw."+year+".nc","r")
    time = cd_calendar(f->time, 2)
    ind_start = ind(time.eq.start_date_LaNina(i))
    ind_end   = ind(time.eq.end_date_LaNina(i))
    ndays_i   = ind_end - ind_start + 1
    var_LaNina_patdays(idx_pat:idx_pat+ndays_i-1,:,:) = \
        f->$var$(ind_start:ind_end,{level},{lat_s:lat_n},{lon_w:lon_e})
    var_LaNina(i,:,:) = \
        dim_sum_n_Wrap(f->$var$(ind_start:ind_end,{level},{lat_s:lat_n},{lon_w:lon_e}),0)
    idx_pat = idx_pat + ndays_i
    delete([/time, ind_start, ind_end, ndays_i/])
end do
var_LaNina_total        = dim_sum_n_Wrap(var_LaNina,0)
var_LaNina_pattern_mean = var_LaNina_total/LaNina_event_total_days
copy_VarCoords(var_LaNina_total, var_LaNina_pattern_mean)
var_LaNina_pattern_prob = ttest_1sample_lag1(var_LaNina_patdays)
; 清理La Nina中间变量
delete([/var_LaNina, var_LaNina_patdays, var_LaNina_total/])

; ===========================================================
; -----输出 NetCDF
; ===========================================================
out_file = "pattern_out_"+var+".nc"

setfileoption("nc","Format","LargeFile")

if (fileexists(out_file)) then
    system("rm -f "+out_file)
end if
f_out = addfile(out_file,"c")

f_out->var_EPElNino_pattern_mean = var_EPElNino_pattern_mean
f_out->var_EPElNino_pattern_prob = var_EPElNino_pattern_prob
f_out->var_CPElNino_pattern_mean = var_CPElNino_pattern_mean
f_out->var_CPElNino_pattern_prob = var_CPElNino_pattern_prob
f_out->var_LaNina_pattern_mean   = var_LaNina_pattern_mean
f_out->var_LaNina_pattern_prob   = var_LaNina_pattern_prob