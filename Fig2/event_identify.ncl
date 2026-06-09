; ;-----SON define
; EP_ElNino_year = (/ 1983, 1988, 1998, 2007, 2016, 2024/)
; CP_ElNino_year = (/ 1987, 1992, 1995, 2003, 2005, 2010, 2019/)
; ; EP_LaNina_year = (/ 1985, 1996, 2000, 2008, 2018, 2021/)
; LaNina_year = (/ 1984, 1989, 1999, 2001, 2011, 2012, 2022, 2023/)

EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016/)
CP_ElNino_year = (/1987, 1992, 1995, 2003, 2010, 2019/)
EP_LaNina_year = (/1996, 2000, 2008, 2021/)
CP_LaNina_year = (/1984, 1989, 1999, 2011, 2012, 2022, 2023/)
LaNina_year = array_append_record(EP_LaNina_year,CP_LaNina_year,0)

;-----standardrize the entire series SON-SON-...-SON
;-----日期位置索引相对于8.17-12.15！
yearStart = 1980
yearEnd = 2024
years = ispan(yearStart,yearEnd,1)
start_month = 8
start_day = 17
day = ispan(1,121,1) ;17 Aug - SON - 15 Dec
nyear = dimsizes(years)
nday = dimsizes(day)

f = addfile("./pattern_index.nc","r")
EPElNino = f->p_cor_EPElNino
EPElNino_index = onedtond(dim_standardize(ndtooned(EPElNino), 0), (/nyear,nday/))
CPElNino = f->p_cor_CPElNino
CPElNino_index = onedtond(dim_standardize(ndtooned(CPElNino), 0), (/nyear,nday/))
LaNina = f->p_cor_LaNina
LaNina_index = onedtond(dim_standardize(ndtooned(LaNina), 0), (/nyear,nday/))

copy_VarCoords(EPElNino, EPElNino_index)
copy_VarCoords(CPElNino, CPElNino_index)
copy_VarCoords(LaNina, LaNina_index)

threshold = 1.  ; Threshold for extreme event
min_duration = 10  ; Minimum consecutive days

year_out = new((/100/),integer)
year_out = 999
start_idx_out = year_out
end_idx_out = year_out
max_idx_out = year_out

;-----识别主体：持续超过阈值n天的事件，记录起始、峰值和结束索引
m = 0
do year = yearStart,yearEnd
    ;El Nino
    if any(EP_ElNino_year.eq.year).or.any(CP_ElNino_year.eq.year)then
        yr = year-1
        if any(EP_ElNino_year.eq.year) then
            index = EPElNino_index({yr}, :)
        elseif any(CP_ElNino_year.eq.year)
            index = CPElNino_index({yr}, :)
        end if
        ; --- Initialize date array ---
        start_year = yr
        time_units = "days since "+start_year+"-01-01 00:00:00"
        dates = new(nday, integer)
        time = cd_inv_calendar(start_year, start_month, start_day, 0, 0, 0, time_units, 0)
        do i = 0, nday-1
            time_tmp = time + i
            time_tmp@units = time_units
            date = tointeger(ut_calendar(time_tmp, 0))
            dates(i) = date(0,0)*10000+date(0,1)*100+date(0,2)
        end do

        ; --- Detect extreme events ---
        events = new(100, integer, -999)  ; Store start indices of events
        event_count = 0
        i = 0
        do while (i .lt. nday)
            if (index(i) .gt. threshold) then
                ; Start of a potential event
                start_idx = i
                count = 1
                do while (i+1 .lt. nday .and. index(i+1) .gt. threshold)
                    i = i + 1
                    count = count + 1
                end do
                ; Check if duration meets minimum requirement
                if (count .ge. min_duration) then
                    events(event_count) = start_idx
                    event_count = event_count + 1
                end if
                i = i + 1
            else
                i = i + 1
            end if
        end do

        ; --- Step 4: Process and output results ---
        if (event_count .eq. 0) then
            print(yr+": No extreme events found.")
        else
            print(yr+": Found " + event_count + " extreme events:")
            do j = 0, event_count-1
            start_idx = events(j)
            ; Find end of event
            end_idx = start_idx
            do while (end_idx+1 .lt. nday .and. index(end_idx+1) .gt. threshold)
                end_idx = end_idx + 1
            end do
            ; Find maximum index value and its position
            event_indices = ind(start_idx .le. ispan(0, nday-1, 1) .and. ispan(0, nday-1, 1) .le. end_idx)
            max_val = max(index(event_indices))
            max_idx = event_indices(ind(index(event_indices) .eq. max_val))
            max_date = dates(max_idx)
            delete(event_indices)
            year_out(m) = yr
            start_idx_out(m) = start_idx
            end_idx_out(m) = end_idx
            max_idx_out(m) = max_idx
            m = m+1
            end do
        end if
    end if

    ;La Nina
    if any(LaNina_year.eq.year) then
        yr = year-1
        index = LaNina_index({yr}, :)

        ; --- Initialize date array ---
        start_year = yr
        time_units = "days since "+start_year+"-01-01 00:00:00"
        dates = new(nday, integer)
        time = cd_inv_calendar(start_year, start_month, start_day, 0, 0, 0, time_units, 0)
        do i = 0, nday-1
            time_tmp = time + i
            time_tmp@units = time_units
            date = tointeger(ut_calendar(time_tmp, 0))
            dates(i) = date(0,0)*10000+date(0,1)*100+date(0,2)
        end do

        ; --- Detect extreme events ---
        threshold = 1.0  ; Threshold for extreme event
        events = new(100, integer, -999)  ; Store start indices of events
        event_count = 0
        i = 0
        do while (i .lt. nday)
            if (index(i) .gt. threshold) then
                ; Start of a potential event
                start_idx = i
                count = 1
                do while (i+1 .lt. nday .and. index(i+1) .gt. threshold)
                    i = i + 1
                    count = count + 1
                end do
                ; Check if duration meets minimum requirement
                if (count .ge. min_duration) then
                    events(event_count) = start_idx
                    event_count = event_count + 1
                end if
                i = i + 1
            else
                i = i + 1
            end if
        end do

        ; --- Step 4: Process and output results ---
        if (event_count .eq. 0) then
            print(yr+": No extreme events found.")
        else
            print(yr+": Found " + event_count + " extreme events:")
            do j = 0, event_count-1
            start_idx = events(j)
            ; Find end of event
            end_idx = start_idx
            do while (end_idx+1 .lt. nday .and. index(end_idx+1) .gt. threshold)
                end_idx = end_idx + 1
            end do
            ; Find maximum index value and its position
            event_indices = ind(start_idx .le. ispan(0, nday-1, 1) .and. ispan(0, nday-1, 1) .le. end_idx)
            max_val = max(index(event_indices))
            max_idx = event_indices(ind(index(event_indices) .eq. max_val))
            max_date = dates(max_idx)
            delete(event_indices)
            year_out(m) = yr
            start_idx_out(m) = start_idx
            end_idx_out(m) = end_idx
            max_idx_out(m) = max_idx
            m = m+1
            end do
        end if
    end if
end do

;-----至此得到事件年份、起始、峰值和结束索引，3类结果统一存储在一起（便于绘制event_statistic分布图）
;-----去除多余预分配空间
year_out_tofile = year_out(ind(year_out.ne.999))
start_idx_out_tofile = start_idx_out(ind(start_idx_out.ne.999))
end_idx_out_tofile = end_idx_out(ind(end_idx_out.ne.999))
max_idx_out_tofile = max_idx_out(ind(max_idx_out.ne.999))

;-----接下来把结果分类存储，便于后续分析
;-----构建8.17-12.15的时间序列，提取起始、峰值和结束对应的时间
;-----先分别统计相应型的个数EP CP LaNina
nEPElNino = 0
do i = 0, dimsizes(year_out_tofile)-1
    if any(EP_ElNino_year .eq.(year_out_tofile(i)+1)) then
        nEPElNino = nEPElNino + 1
    end if
end do

nCPElNino = 0
do i = 0, dimsizes(year_out_tofile)-1
    if any(CP_ElNino_year .eq.(year_out_tofile(i)+1)) then
        nCPElNino = nCPElNino + 1
    end if
end do

nLaNina = 0
do i = 0, dimsizes(year_out_tofile)-1
    if any(LaNina_year .eq.(year_out_tofile(i)+1)) then
        nLaNina = nLaNina + 1
    end if
end do
print("nEPElNino: "+nEPElNino)
print("nCPElNino: "+nCPElNino)
print("nLaNina: "+nLaNina)

start_date_EPElNino = new(nEPElNino, integer)
end_date_EPElNino = new(nEPElNino, integer)
max_date_EPElNino = new(nEPElNino, integer)
start_date_CPElNino = new(nCPElNino, integer)
end_date_CPElNino = new(nCPElNino, integer)
max_date_CPElNino = new(nCPElNino, integer)
start_date_LaNina = new(nLaNina, integer)
end_date_LaNina = new(nLaNina, integer)
max_date_LaNina = new(nLaNina, integer)

start_idx_EPElNino = new(nEPElNino, integer)
end_idx_EPElNino = new(nEPElNino, integer)
max_idx_EPElNino = new(nEPElNino, integer)
start_idx_CPElNino = new(nCPElNino, integer)
end_idx_CPElNino = new(nCPElNino, integer)
max_idx_CPElNino = new(nCPElNino, integer)
start_idx_LaNina = new(nLaNina, integer)
end_idx_LaNina = new(nLaNina, integer)
max_idx_LaNina = new(nLaNina, integer)


j = 0 ;EP CP La循环指针
m = 0
n = 0
do i = 0,dimsizes(year_out_tofile)-1
    yr = year_out_tofile(i)
    ;逐年建立8.17-12.15时间序列
    start_year = yr
    time_units = "days since "+start_year+"-01-01 00:00:00"
    dates = new(nday, integer)
    time = cd_inv_calendar(start_year, start_month, start_day, 0, 0, 0, time_units, 0)
    do k = 0, nday-1
        time_tmp = time + k
        time_tmp@units = time_units
        date = tointeger(ut_calendar(time_tmp, 0))
        dates(k) = date(0,0)*10000+date(0,1)*100+date(0,2)
    end do

    ;根据年份归类
    if any(EP_ElNino_year .eq.(yr+1)) then
        start_date_EPElNino(j) = dates(start_idx_out_tofile(i))
        end_date_EPElNino(j) = dates(end_idx_out_tofile(i))
        max_date_EPElNino(j) = dates(max_idx_out_tofile(i))
        start_idx_EPElNino(j) = start_idx_out_tofile(i)
        end_idx_EPElNino(j) = end_idx_out_tofile(i)
        max_idx_EPElNino(j) = max_idx_out_tofile(i)
        j = j + 1
    elseif any(CP_ElNino_year .eq.(yr+1)) then
        start_date_CPElNino(m) = dates(start_idx_out_tofile(i))
        end_date_CPElNino(m) = dates(end_idx_out_tofile(i))
        max_date_CPElNino(m) = dates(max_idx_out_tofile(i))
        start_idx_CPElNino(m) = start_idx_out_tofile(i)
        end_idx_CPElNino(m) = end_idx_out_tofile(i)
        max_idx_CPElNino(m) = max_idx_out_tofile(i)
        m = m + 1
    elseif any(LaNina_year .eq.(yr+1)) then
        start_date_LaNina(n) = dates(start_idx_out_tofile(i))
        end_date_LaNina(n) = dates(end_idx_out_tofile(i))
        max_date_LaNina(n) = dates(max_idx_out_tofile(i))
        start_idx_LaNina(n) = start_idx_out_tofile(i)
        end_idx_LaNina(n) = end_idx_out_tofile(i)
        max_idx_LaNina(n) = max_idx_out_tofile(i)
        n = n + 1
    end if
end do

;-----output the Comp_avg_hgt
	out_file = "event_identify.nc"
	if (fileexists(out_file)) then
		system("rm -f "+out_file)
	end if
	f_out = addfile(out_file,"c")
    f_out->year = year_out_tofile
    f_out->start_idx = start_idx_out_tofile
    f_out->end_idx = end_idx_out_tofile
	f_out->max_idx = max_idx_out_tofile
    f_out->nEPElNino = nEPElNino
    f_out->nCPElNino = nCPElNino
    f_out->nLaNina = nLaNina
    f_out->start_date_EPElNino = start_date_EPElNino
    f_out->end_date_EPElNino = end_date_EPElNino
    f_out->max_date_EPElNino = max_date_EPElNino
    f_out->start_date_CPElNino = start_date_CPElNino
    f_out->end_date_CPElNino = end_date_CPElNino
    f_out->max_date_CPElNino = max_date_CPElNino
    f_out->start_date_LaNina = start_date_LaNina
    f_out->end_date_LaNina = end_date_LaNina
    f_out->max_date_LaNina = max_date_LaNina
    f_out->start_idx_EPElNino = start_idx_EPElNino
    f_out->end_idx_EPElNino = end_idx_EPElNino
    f_out->max_idx_EPElNino = max_idx_EPElNino
    f_out->start_idx_CPElNino = start_idx_CPElNino
    f_out->end_idx_CPElNino = end_idx_CPElNino
    f_out->max_idx_CPElNino = max_idx_CPElNino
    f_out->start_idx_LaNina = start_idx_LaNina
    f_out->end_idx_LaNina = end_idx_LaNina
    f_out->max_idx_LaNina = max_idx_LaNina