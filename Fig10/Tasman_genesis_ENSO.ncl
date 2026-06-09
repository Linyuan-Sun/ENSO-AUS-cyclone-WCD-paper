
;==================== 用户设置部分 ====================
dir       = "./"  ; 输入文件路径
in_file = "Tasman_genesis.csv"
out_file1  = "Tasman_genesis_EPElNino.csv"             ; 输出文件名
out_file2  = "Tasman_genesis_CPElNino.csv"             ; 输出文件名
out_file3  = "Tasman_genesis_LaNina.csv"               ; 输出
;======================================================

;---加载型事件日期
f = addfile("/g/data/m35/ls3248/ENSO_cyclone/paper/Fig2/event_identify.nc", "r")
start_date_EPElNino = f->start_date_EPElNino
end_date_EPElNino = f->end_date_EPElNino
start_date_CPElNino = f->start_date_CPElNino
end_date_CPElNino = f->end_date_CPElNino
start_date_LaNina = f->start_date_LaNina
end_date_LaNina = f->end_date_LaNina

;---单独处理EP ElNino缺失2023年
EP_missing = 2024-1
EPElNino_event_year = start_date_EPElNino/10000
start_date_EPElNino := start_date_EPElNino(ind(EPElNino_event_year.ne.EP_missing))
end_date_EPElNino := end_date_EPElNino(ind(EPElNino_event_year.ne.EP_missing))

;--- 若输出文件存在则删除
if (fileexists(out_file1)) then
    system("rm -f " + out_file1)
end if
if (fileexists(out_file2)) then
    system("rm -f " + out_file2)
end if
if (fileexists(out_file3)) then
    system("rm -f " + out_file3)
end if

;--- 写入表头
system("echo " + "datetime,lon,lat,pmin,label" + " >> " + out_file1)
system("echo " + "datetime,lon,lat,pmin,label" + " >> " + out_file2)
system("echo " + "datetime,lon,lat,pmin,label" + " >> " + out_file3)

;---读入Tasman轨迹文件
fname = dir+in_file                              
lines = asciiread(fname, -1, "string")                           ; 读取 CSV 文件的所有行到字符串数组
track_info = lines(1:)                                           ; 跳过表头，提取数据行
datetime_str = str_get_field(track_info,1,",")                   ; 提取第1列的日期时间字段并转换为YYYYMMDD整型时间数组
datetime = toint(str_get_cols(datetime_str, 0, 7))               ; 提取前8位 YYYYMMDD
label    = tointeger(str_get_field(track_info, 5, ","))

;--- 唯一气旋编号
uniq_label = get_unique_values(label)
nlab = dimsizes(uniq_label)

;--- 遍历每个气旋
do ilab = 0, nlab-1
    this_label = uniq_label(ilab)
    idx = ind(label.eq.this_label)
    track_info_cyclone = track_info(idx)
    datetime_cyclone = datetime(idx)
    first_datetime = datetime_cyclone(0);--- 获取该气旋的生成日期

    ; 是否属于EPElNino事件期
    do i = 0, dimsizes(start_date_EPElNino)-1
        if (first_datetime.ge.start_date_EPElNino(i) .and. first_datetime.le.end_date_EPElNino(i)) then
            write_table(out_file1, "a", [/track_info_cyclone/], "%s")    ; 将符合条件的轨迹信息追加写入输出文件
        end if
    end do

    ; 是否属于CPElNino事件期
    do i = 0, dimsizes(start_date_CPElNino)-1
        if (first_datetime.ge.start_date_CPElNino(i) .and. first_datetime.le.end_date_CPElNino(i)) then
            write_table(out_file2, "a", [/track_info_cyclone/], "%s")    ; 将符合条件的轨迹信息追加写入输出文件
        end if
    end do

    ; 是否属于LaNina事件期
    do i = 0, dimsizes(start_date_LaNina)-1
        if (first_datetime.ge.start_date_LaNina(i) .and. first_datetime.le.end_date_LaNina(i)) then
            write_table(out_file3, "a", [/track_info_cyclone/], "%s")    ; 将符合条件的轨迹信息追加写入输出文件
        end if
    end do

    delete([/idx,track_info_cyclone,datetime_cyclone/])
end do
