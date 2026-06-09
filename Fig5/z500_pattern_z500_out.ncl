pattern = "z500"
var = "z"
level = 500
nlat = 180
nlon = 360
season = "SON"

delta_day = 10

; load event index
f = addfile("/g/data/m35/ls3248/ENSO_synop/comp_pattern_based/event_iden/event_identify.nc", "r")
event_year_temp = f->year
start_idx = f->start_idx

; 去掉缺省值999
event_year = event_year_temp(ind(event_year_temp.ne.999))
event_year_num = dimsizes(event_year)

EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016, 2024/)
CP_ElNino_year = (/1987, 1992, 1995, 2003, 2005, 2010, 2019/)
EP_LaNina_year = (/1996, 2000, 2008, 2021/)
CP_LaNina_year = (/1984, 1989, 1999, 2011, 2012, 2022, 2023/)
LaNina_year = array_append_record(EP_LaNina_year,CP_LaNina_year,0)

; 获取相应类型，先计数，再赋值
; EP El Nino
k = 0
do i = 0,event_year_num-1
    if any(EP_ElNino_year-1.eq.event_year(i))then
        k = k+1
    end if
end do
event_year_EPElNino = new((/k/),integer)
start_idx_EPElNino = new((/k/),integer)
m = 0
do i = 0,event_year_num-1
    if any(EP_ElNino_year-1.eq.event_year(i))then
        event_year_EPElNino(m) = event_year(i)
        start_idx_EPElNino(m) = start_idx(i)
        m = m+1
    end if
end do

; CP El Nino
k = 0
do i = 0,event_year_num-1
    if any(CP_ElNino_year-1.eq.event_year(i))then
        k = k+1
    end if
end do
event_year_CPElNino = new((/k/),integer)
start_idx_CPElNino = new((/k/),integer)
m = 0
do i = 0,event_year_num-1
    if any(CP_ElNino_year-1.eq.event_year(i))then
        event_year_CPElNino(m) = event_year(i)
        start_idx_CPElNino(m) = start_idx(i)
        m = m+1
    end if
end do

; La Nina
k = 0
do i = 0,event_year_num-1
    if any(LaNina_year-1.eq.event_year(i))then
        k = k+1
    end if
end do
event_year_LaNina = new((/k/),integer)
start_idx_LaNina = new((/k/),integer)
m = 0
do i = 0,event_year_num-1
    if any(LaNina_year-1.eq.event_year(i))then
        event_year_LaNina(m) = event_year(i)
        start_idx_LaNina(m) = start_idx(i)
        m = m+1
    end if
end do
; printVarSummary(event_year_EPElNino)
; printVarSummary(start_idx_EPElNino)
; printVarSummary(event_year_CPElNino)
; printVarSummary(start_idx_CPElNino)
; printVarSummary(event_year_LaNina)
; printVarSummary(start_idx_LaNina)

;--------------------------
;------load var data-------
;--------------------------
;EP El Nino
k = dimsizes(event_year_EPElNino)
var_EPElNino = new((/k,delta_day*2+1,nlat,nlon/),"float") ;-delta_day max +delta_day
do i = 0,k-1
    if start_idx_EPElNino(i).ge.0+delta_day .and. start_idx_EPElNino(i).le.120-delta_day then ;为保证最值前后有+-delta_day天的数据，最值在窗口以外的事件不计入
        year = event_year_EPElNino(i)
        f = addfile("/g/data/if69/ls3248/data/ERA5/daily/"+var+"/LF/"+var+".ERA5."+season+".anom.LF."+year+".nc","r")
        var_EPElNino(i,:,:,:) = f->$var$(start_idx_EPElNino(i)-delta_day:start_idx_EPElNino(i)+delta_day,{level},:,:)
    end if
end do
print("EP ElNino is done.")

;CP El Nino
k = dimsizes(event_year_CPElNino)
var_CPElNino = new((/k,delta_day*2+1,nlat,nlon/),"float") ;-delta_day max +delta_day
do i = 0,k-1
    if start_idx_CPElNino(i).ge.0+delta_day .and. start_idx_CPElNino(i).le.120-delta_day then ;为保证最值前后有+-delta_day天的数据，最值在窗口以外的事件不计入
        year = event_year_CPElNino(i)
        f = addfile("/g/data/if69/ls3248/data/ERA5/daily/"+var+"/LF/"+var+".ERA5."+season+".anom.LF."+year+".nc","r")
        var_CPElNino(i,:,:,:) = f->$var$(start_idx_CPElNino(i)-delta_day:start_idx_CPElNino(i)+delta_day,{level},:,:)
    end if
end do
print("CP ElNino is done.")

;La Nina
k = dimsizes(event_year_LaNina)
var_LaNina = new((/k,delta_day*2+1,nlat,nlon/),"float") ;-delta_day max +delta_day
do i = 0,k-1
    if start_idx_LaNina(i).ge.0+delta_day .and. start_idx_LaNina(i).le.120-delta_day then ;为保证最值前后有+-delta_day天的数据，最值在窗口以外的事件不计入
        year = event_year_LaNina(i)
        f = addfile("/g/data/if69/ls3248/data/ERA5/daily/"+var+"/LF/"+var+".ERA5."+season+".anom.LF."+year+".nc","r")
        var_LaNina(i,:,:,:) = f->$var$(start_idx_LaNina(i)-delta_day:start_idx_LaNina(i)+delta_day,{level},:,:)
    end if
end do
print("LaNina is done.")

;-----output the Comp_avg_hgt
out_file = pattern+"_pattern_"+var+level+"_"+season+".nc"

varout1 = var+"_EPElNino"
varout3 = var+"_CPElNino"
varout4 = var+"_LaNina"

setfileoption("nc","Format","LargeFile")

if (fileexists(out_file)) then
    system("rm -f "+out_file)
end if
f_out = addfile(out_file,"c")
f_out->$varout1$ = var_EPElNino
f_out->$varout3$ = var_CPElNino
f_out->$varout4$ = var_LaNina