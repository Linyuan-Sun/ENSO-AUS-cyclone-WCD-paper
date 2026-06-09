
;==================== 用户设置部分 ====================
dir       = "/g/data/if69/ls3248/data/weather/cyclone/track/"  ; 输入文件路径
out_file  = "Tasman_genesis.csv"             ; 输出文件名
lon_min   = 153.0                            ; 塔斯曼海经度范围
lon_max   = 166.0
lat_min   = -45.0                            ; 纬度为负值（南纬）
lat_max   = -30.0
years     = ispan(1979, 2022, 1)             ; 1979–2022 年
months    = (/8,9,10,11,12/)                 ; ASOND
start_md = 817   ; 8月17日
end_md = 1215    ; 12月15日
min_lifetime_points = 12  ; 最小生命期点数（每小时1个点，12表示12小时）
;======================================================

;--- 若输出文件存在则删除
if (fileexists(out_file)) then
system("rm -f " + out_file)
end if

;--- 写入表头
system("echo " + "datetime,lon,lat,pmin,label" + " >> " + out_file)

;--- 外层循环：年份
do yy = 0, dimsizes(years)-1
year = years(yy)

    ;--- 内层循环：月份
    do mm = 0, dimsizes(months)-1
        month = months(mm)
        fname = dir + "tr_" + sprinti("%04i", year) + sprinti("%02i", month) + ".csv"

        print(">>> 正在处理: " + fname)

        ;--- 读取文件为字符串
        lines = asciiread(fname, -1, "string")
        track_info = lines(1:)

        ;--- 使用str_get_field解析各列
        datetime_str = str_get_field(track_info, 1, ",")
        datetime = toint(str_get_cols(datetime_str, 0, 7))     
        lon      = tofloat(str_get_field(track_info, 2, ","))
        lat      = tofloat(str_get_field(track_info, 3, ","))
        pmin     = tofloat(str_get_field(track_info, 4, ","))
        label    = tointeger(str_get_field(track_info, 5, ","))

        ;--- 唯一气旋编号
        uniq_label = get_unique_values(label)
        nlab = dimsizes(uniq_label)

        ;--- 遍历每个气旋
        do ilab = 0, nlab-1
            this_label = uniq_label(ilab)
            idx = ind(label.eq.this_label)
            ;--- 获取该气旋的生成点（第一个点）
            first_lon = lon(idx(0))
            first_lat = lat(idx(0))
            first_datetime = datetime(idx(0))

            npts = dimsizes(idx)                     ; 该气旋的总记录点数（每小时1个）

            ;========== 新增：生命期筛选 ==========
            if (npts .lt. min_lifetime_points) then
                delete(idx)
                continue                             ; 少于一定时间直接跳过
            end if
            ;======================================

            ;=== 判断生成点是否在塔斯曼海区域且限定时间窗口内 ===
            season_start = year*10000+start_md 
            season_end   = year*10000+end_md
            if (first_lon.ge.lon_min .and. first_lon.le.lon_max .and. first_lat.ge.lat_min .and. first_lat.le.lat_max .and. first_datetime.ge.season_start .and. first_datetime.le.season_end) then
                track_out = track_info(idx)
                write_table(out_file, "a", [/track_out/], "%s")    ; 将符合条件的轨迹信息追加写入输出文件
                delete(track_out)
            end if
            delete(idx)
        end do
        delete([/lines, track_info, datetime_str, datetime, lon, lat, pmin, label, uniq_label/])
    end do  ; month
end do      ; year

print("所有年份和月份处理完成，结果输出至: " + out_file)