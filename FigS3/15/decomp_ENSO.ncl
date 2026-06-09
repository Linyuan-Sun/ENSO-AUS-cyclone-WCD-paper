;-----提取pattern days信息
; ENSO events
EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016/);, 2024/)
EP_missing = 2023

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

;-----单独处理EP ElNino缺失2023年
EPElNino_event_year = start_date_EPElNino/10000
nEPElNino = nEPElNino - num(EPElNino_event_year.eq.EP_missing)
start_date_EPElNino := start_date_EPElNino(ind(EPElNino_event_year.ne.EP_missing))
end_date_EPElNino := end_date_EPElNino(ind(EPElNino_event_year.ne.EP_missing))
start_idx_EPElNino := start_idx_EPElNino(ind(EPElNino_event_year.ne.EP_missing))
end_idx_EPElNino := end_idx_EPElNino(ind(EPElNino_event_year.ne.EP_missing))

;-----计算各事件窗口的总天数
EPElNino_event_total_days = sum(end_idx_EPElNino - start_idx_EPElNino + 1) 
CPElNino_event_total_days = sum(end_idx_CPElNino - start_idx_CPElNino + 1)
LaNina_event_total_days = sum(end_idx_LaNina - start_idx_LaNina + 1)

;气旋日信息
f = addfile("cyclone_flag.nc","r")
flag_temp = f->flag

;AGCD降水维度
nlat = 70
nlon = 89

;-----EPElNino
print("EPElNino:-----")
;型事件期间内的变量和
cyclone_precip_EPElNino = new((/nlat,nlon/),"float")
cyclone_precip_EPElNino = 0.
cyclone_days_num_EPElNino = 0.
noncyclone_precip_EPElNino = new((/nlat,nlon/),"float")
noncyclone_precip_EPElNino = 0.
noncyclone_days_num_EPElNino = 0.
do i = 0,nEPElNino-1
    year = start_date_EPElNino(i)/10000
    start_date = start_date_EPElNino(i)-year*10000
    end_date = end_date_EPElNino(i)-year*10000
    ;-----气旋日
    cyclone_ind = ind(flag_temp({year},{start_date:end_date}).eq.1) ;确定事件内的气旋日索引
    cyclone_num = num(flag_temp({year},{start_date:end_date}).eq.1)
    if cyclone_num .eq. 0 then
        print("No cyclone day during "+start_date_EPElNino(i)+" "+end_date_EPElNino(i))
    else 
        flag_temp_temp = flag_temp({year},{start_date:end_date})
        days = flag_temp_temp&day
        cyclone_date = year*10000+days(cyclone_ind)
        cyclone_days_num_EPElNino = cyclone_days_num_EPElNino + cyclone_num
        ; 读取AGCD降水数据
        f = addfile("/g/data/if69/ls3248/data/AGCD/precip/daily/agcd_v1_precip_total_r05_daily_"+year+".nc","r")
        precip = f->precip
        precip!0 = "time"
        precip&time = cd_calendar(f->time, 2)
        ;识别出1天或其他
        if dimsizes(cyclone_date) .eq. 1 then
            cyclone_precip_EPElNino = cyclone_precip_EPElNino + precip({cyclone_date},:,:)
        else
            cyclone_precip_EPElNino = cyclone_precip_EPElNino + dim_sum_n_Wrap(precip({cyclone_date},:,:),0)
        end if
        delete([/days,flag_temp_temp,precip,cyclone_date/])
    end if
    ;-----非气旋日
    noncyclone_ind = ind(flag_temp({year},{start_date:end_date}).eq.0) ;确定事件内的非气旋日索引
    noncyclone_num = num(flag_temp({year},{start_date:end_date}).eq.0)
    if noncyclone_num .eq. 0 then
        print("No non-cyclone day during "+start_date_EPElNino(i)+" "+end_date_EPElNino(i))
    else 
        flag_temp_temp = flag_temp({year},{start_date:end_date})
        days = flag_temp_temp&day
        noncyclone_date = year*10000+days(noncyclone_ind)
        noncyclone_days_num_EPElNino = noncyclone_days_num_EPElNino + noncyclone_num
        ; 读取AGCD降水数据
        f = addfile("/g/data/if69/ls3248/data/AGCD/precip/daily/agcd_v1_precip_total_r05_daily_"+year+".nc","r")
        precip = f->precip
        precip!0 = "time"
        precip&time = cd_calendar(f->time, 2)
        ;识别出1天或其他
        if dimsizes(noncyclone_date) .eq. 1 then
            noncyclone_precip_EPElNino = noncyclone_precip_EPElNino + precip({noncyclone_date},:,:)
        else
            noncyclone_precip_EPElNino = noncyclone_precip_EPElNino + dim_sum_n_Wrap(precip({noncyclone_date},:,:),0)
        end if
        delete([/days,flag_temp_temp,precip,noncyclone_date/])
    end if
    delete([/cyclone_ind,noncyclone_ind/])
    print(start_date_EPElNino(i)+"-"+end_date_EPElNino(i)+": DONE.")
end do
pc_EPElNino = cyclone_precip_EPElNino/(cyclone_days_num_EPElNino+noncyclone_days_num_EPElNino)
pnc_EPElNino = noncyclone_precip_EPElNino/(cyclone_days_num_EPElNino+noncyclone_days_num_EPElNino)

;-----CPElNino
print("CPElNino:-----")
;型事件期间内的变量和
cyclone_precip_CPElNino = new((/nlat,nlon/),"float")
cyclone_precip_CPElNino = 0.
cyclone_days_num_CPElNino = 0.
noncyclone_precip_CPElNino = new((/nlat,nlon/),"float")
noncyclone_precip_CPElNino = 0.
noncyclone_days_num_CPElNino = 0.
do i = 0,nCPElNino-1
    year = start_date_CPElNino(i)/10000
    start_date = start_date_CPElNino(i)-year*10000
    end_date = end_date_CPElNino(i)-year*10000
    ;-----气旋日
    cyclone_ind = ind(flag_temp({year},{start_date:end_date}).eq.1) ;确定事件内的气旋日索引
    cyclone_num = num(flag_temp({year},{start_date:end_date}).eq.1)
    if cyclone_num .eq. 0 then
        print("No cyclone day during "+start_date_CPElNino(i)+" "+end_date_CPElNino(i))
    else 
        flag_temp_temp = flag_temp({year},{start_date:end_date})
        days = flag_temp_temp&day
        cyclone_date = year*10000+days(cyclone_ind)
        cyclone_days_num_CPElNino = cyclone_days_num_CPElNino + cyclone_num
        ; 读取AGCD降水数据
        f = addfile("/g/data/if69/ls3248/data/AGCD/precip/daily/agcd_v1_precip_total_r05_daily_"+year+".nc","r")
        precip = f->precip
        precip!0 = "time"
        precip&time = cd_calendar(f->time, 2)
        ;识别出1天或其他
        if dimsizes(cyclone_date) .eq. 1 then
            cyclone_precip_CPElNino = cyclone_precip_CPElNino + precip({cyclone_date},:,:)
        else
            cyclone_precip_CPElNino = cyclone_precip_CPElNino + dim_sum_n_Wrap(precip({cyclone_date},:,:),0)
        end if
        delete([/days,flag_temp_temp,precip,cyclone_date/])
    end if
    ;-----非气旋日
    noncyclone_ind = ind(flag_temp({year},{start_date:end_date}).eq.0) ;确定事件内的非气旋日索引
    noncyclone_num = num(flag_temp({year},{start_date:end_date}).eq.0)
    if noncyclone_num .eq. 0 then
        print("No non-cyclone day during "+start_date_CPElNino(i)+" "+end_date_CPElNino(i))
    else 
        flag_temp_temp = flag_temp({year},{start_date:end_date})
        days = flag_temp_temp&day
        noncyclone_date = year*10000+days(noncyclone_ind)
        noncyclone_days_num_CPElNino = noncyclone_days_num_CPElNino + noncyclone_num
        ; 读取AGCD降水数据
        f = addfile("/g/data/if69/ls3248/data/AGCD/precip/daily/agcd_v1_precip_total_r05_daily_"+year+".nc","r")
        precip = f->precip
        precip!0 = "time"
        precip&time = cd_calendar(f->time, 2)
        ;识别出1天或其他
        if dimsizes(noncyclone_date) .eq. 1 then
            noncyclone_precip_CPElNino = noncyclone_precip_CPElNino + precip({noncyclone_date},:,:)
        else
            noncyclone_precip_CPElNino = noncyclone_precip_CPElNino + dim_sum_n_Wrap(precip({noncyclone_date},:,:),0)
        end if
        delete([/days,flag_temp_temp,precip,noncyclone_date/])
    end if
    delete([/cyclone_ind,noncyclone_ind/])
    print(start_date_CPElNino(i)+"-"+end_date_CPElNino(i)+": DONE.")
end do
pc_CPElNino = cyclone_precip_CPElNino/(cyclone_days_num_CPElNino+noncyclone_days_num_CPElNino)
pnc_CPElNino = noncyclone_precip_CPElNino/(cyclone_days_num_CPElNino+noncyclone_days_num_CPElNino)

;-----LaNina
print("LaNina:-----")
;型事件期间内的变量和
cyclone_precip_LaNina = new((/nlat,nlon/),"float")
cyclone_precip_LaNina = 0.
cyclone_days_num_LaNina = 0.
noncyclone_precip_LaNina = new((/nlat,nlon/),"float")
noncyclone_precip_LaNina = 0.
noncyclone_days_num_LaNina = 0.
do i = 0,nLaNina-1
    year = start_date_LaNina(i)/10000
    start_date = start_date_LaNina(i)-year*10000
    end_date = end_date_LaNina(i)-year*10000
    ;-----气旋日
    cyclone_ind = ind(flag_temp({year},{start_date:end_date}).eq.1) ;确定事件内的气旋日索引
    cyclone_num = num(flag_temp({year},{start_date:end_date}).eq.1)
    if cyclone_num .eq. 0 then
        print("No cyclone day during "+start_date_LaNina(i)+" "+end_date_LaNina(i))
    else 
        flag_temp_temp = flag_temp({year},{start_date:end_date})
        days = flag_temp_temp&day
        cyclone_date = year*10000+days(cyclone_ind)
        cyclone_days_num_LaNina = cyclone_days_num_LaNina + cyclone_num
        ; 读取AGCD降水数据
        f = addfile("/g/data/if69/ls3248/data/AGCD/precip/daily/agcd_v1_precip_total_r05_daily_"+year+".nc","r")
        precip = f->precip
        precip!0 = "time"
        precip&time = cd_calendar(f->time, 2)
        ;识别出1天或其他
        if dimsizes(cyclone_date) .eq. 1 then
            cyclone_precip_LaNina = cyclone_precip_LaNina + precip({cyclone_date},:,:)
        else
            cyclone_precip_LaNina = cyclone_precip_LaNina + dim_sum_n_Wrap(precip({cyclone_date},:,:),0)
        end if
        delete([/days,flag_temp_temp,precip,cyclone_date/])
    end if
    ;-----非气旋日
    noncyclone_ind = ind(flag_temp({year},{start_date:end_date}).eq.0) ;确定事件内的非气旋日索引
    noncyclone_num = num(flag_temp({year},{start_date:end_date}).eq.0)
    if noncyclone_num .eq. 0 then
        print("No non-cyclone day during "+start_date_LaNina(i)+" "+end_date_LaNina(i))
    else 
        flag_temp_temp = flag_temp({year},{start_date:end_date})
        days = flag_temp_temp&day
        noncyclone_date = year*10000+days(noncyclone_ind)
        noncyclone_days_num_LaNina = noncyclone_days_num_LaNina + noncyclone_num
        ; 读取AGCD降水数据
        f = addfile("/g/data/if69/ls3248/data/AGCD/precip/daily/agcd_v1_precip_total_r05_daily_"+year+".nc","r")
        precip = f->precip
        precip!0 = "time"
        precip&time = cd_calendar(f->time, 2)
        ;识别出1天或其他
        if dimsizes(noncyclone_date) .eq. 1 then
            noncyclone_precip_LaNina = noncyclone_precip_LaNina + precip({noncyclone_date},:,:)
        else
            noncyclone_precip_LaNina = noncyclone_precip_LaNina + dim_sum_n_Wrap(precip({noncyclone_date},:,:),0)
        end if
        delete([/days,flag_temp_temp,precip,noncyclone_date/])
    end if
    delete([/cyclone_ind,noncyclone_ind/])
    print(start_date_LaNina(i)+"-"+end_date_LaNina(i)+": DONE.")
end do
pc_LaNina = cyclone_precip_LaNina/(cyclone_days_num_LaNina+noncyclone_days_num_LaNina)
pnc_LaNina = noncyclone_precip_LaNina/(cyclone_days_num_LaNina+noncyclone_days_num_LaNina)

;赋予维度
f = addfile("/g/data/if69/ls3248/data/AGCD/precip/daily/agcd_v1_precip_total_r05_daily_2000.nc","r")
precip = f->precip(0,:,:)
copy_VarCoords(precip, pc_EPElNino)
copy_VarCoords(precip, pnc_EPElNino)
copy_VarCoords(precip, pc_CPElNino)
copy_VarCoords(precip, pnc_CPElNino)
copy_VarCoords(precip, pc_LaNina)
copy_VarCoords(precip, pnc_LaNina)

;-----output the Comp_avg_hgt
out_file = "decomp_ENSO.nc"

if (fileexists(out_file)) then
    system("rm -f "+out_file)
end if
f_out = addfile(out_file,"c")
f_out->pc_EPElNino = pc_EPElNino
f_out->pnc_EPElNino = pnc_EPElNino
f_out->cyclone_days_num_EPElNino = cyclone_days_num_EPElNino
f_out->noncyclone_days_num_EPElNino = noncyclone_days_num_EPElNino

f_out->pc_CPElNino = pc_CPElNino
f_out->pnc_CPElNino = pnc_CPElNino
f_out->cyclone_days_num_CPElNino = cyclone_days_num_CPElNino
f_out->noncyclone_days_num_CPElNino = noncyclone_days_num_CPElNino

f_out->pc_LaNina = pc_LaNina
f_out->pnc_LaNina = pnc_LaNina
f_out->cyclone_days_num_LaNina = cyclone_days_num_LaNina
f_out->noncyclone_days_num_LaNina = noncyclone_days_num_LaNina