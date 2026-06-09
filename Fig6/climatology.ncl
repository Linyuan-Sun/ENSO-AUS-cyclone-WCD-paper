
  ;==== 设置路径与参数 ====
  dir = "/g/data/if69/ls3248/data/weather/cyclone/track/"           ; 设置存放 CSV 文件的目录
  years = ispan(1979,2022,1)                                            ; 创建年份数组，从1979到2022，步长为1
  season = "SON"                                                       ; 定义季节为秋季（SON）
  months = (/9,10,11/)                                                  ; 定义秋季（SON）月份：9月、10月、11月
  nmonth = (/30,31,30/)
  lonMin = -180                                                         ; 经度范围最小值：-180°
  lonMax = 180                                                          ; 经度范围最大值：180°
  latMin = -90                                                          ; 纬度范围最小值：-90°
  latMax = 0                                                            ; 纬度范围最大值：0°（南半球）
  dlon = 5                                                            ; 经度网格间距：2.5°
  dlat = 5                                                           ; 纬度网格间距：2.5°
  lon_bins = fspan(lonMin,lonMax, toint(((lonMax-lonMin)/dlon)+1))     ; 创建经度网格边界数组，[-180, -177.5, ..., 180]
  lat_bins = fspan(latMin,latMax, toint(((latMax-latMin)/dlat)+1))     ; 创建纬度网格边界数组，[-90, -87.5, ..., 0]
  nlon = dimsizes(lon_bins)-1                                           ; 计算经度网格数量（边界数-1）
  nlat = dimsizes(lat_bins)-1                                           ; 计算纬度网格数量（边界数-1）
  track_count = new((/nlat,nlon/), float)                                ; 创建二维数组（纬度×经度）用于记录轨迹点计数
  track_count = 0.0                                                     ; 初始化轨迹计数数组为0

  ;==== 计算总小时数 ====
  total_hours = dimsizes(years) * sum(nmonth) * 24                                ; 总小时数 = 44年 × (9月30天 + 10月31天 + 11月30天) × 24小时 = 96,096小时

  ;==== 循环读取文件 ====
  do yy = 0, dimsizes(years)-1                                          ; 遍历年份数组（1979到2022）
    year = years(yy)                                                    ; 获取当前年份
    do mm = 0, dimsizes(months)-1                                     ; 遍历月份数组（9、10、11月）
      month = months(mm)                                              ; 获取当前月份
      fname = dir + "tr_" + sprinti("%04i", year) + sprinti("%02i", month) + ".csv"  ; 构造文件名，如 tr_197909.csv
      print("Processing file: " + fname)                             ; 打印当前处理的文件名
      if (fileexists(fname)) then                                    ; 检查文件是否存在
        lines = asciiread(fname, -1, "string")                       ; 读取 CSV 文件的所有行到字符串数组
        nline = dimsizes(lines)                                      ; 获取文件行数
        if (nline .gt. 1) then                                      ; 确保文件至少有2行（跳过标题行）
          lon = tofloat(str_get_field(lines(1:),2,","))           ; 提取第2列（经度），从第2行开始，转换为浮点数
          lat = tofloat(str_get_field(lines(1:),3,","))           ; 提取第3列（纬度），从第2行开始，转换为浮点数
          do n = 0, dimsizes(lon)-1                                 ; 遍历每个轨迹点
            if (.not.ismissing(lon(n)) .and. .not.ismissing(lat(n))) then  ; 检查经纬度是否非缺失值
              if (lon(n).ge.lonMin .and. lon(n).lt.lonMax .and. lat(n).ge.latMin .and. lat(n).lt.latMax) then  ; 检查经纬度是否在范围内
                ilon = ind(lon_bins.le.lon(n))                       ; 找到满足 lon_bins <= lon(n) 的索引
                ilat = ind(lat_bins.le.lat(n))                       ; 找到满足 lat_bins <= lat(n) 的索引
                if (dimsizes(ilon).gt.0 .and. dimsizes(ilat).gt.0) then  ; 确保找到有效网格索引
                  track_count(ilat(dimsizes(ilat)-1), ilon(dimsizes(ilon)-1)) = track_count(ilat(dimsizes(ilat)-1), ilon(dimsizes(ilon)-1)) + 1  ; 在对应网格计数+1
                  delete([/ilon,ilat/])
                end if
              end if
            end if
          end do
        end if
      else
        print("Warning: File " + fname + " does not exist, skipping...")  ; 打印文件不存在的警告
      end if
      delete([/lines,lon,lat/])                                               ; 删除读取的行以释放内存
    end do
  end do

    ;==== 转换为频率 ====
    track_density = track_count / tofloat(total_hours)*100                     ; 计算轨迹密度：计数除以总小时数（每小时每网格的频率）
    track_density!0 = "lat"
    track_density!1 = "lon"
    lat_track = lat_bins(0:nlat-1) + dlat/2  ; 纬度网格中心：[-88.75, -86.25, ..., -3.75]
    lon_track = lon_bins(0:nlon-1) + dlon/2  ; 经度网格中心：[-178.75, -176.25, ..., 178.75]
    lat_track@units = "degrees_north"
    lon_track@units = "degrees_east"
    track_density&lat = lat_track
    track_density&lon = lon_track

    ; track_density = where(track_density.gt.0.15,0,track_density)

    ;-----output the Comp_avg_hgt
    out_file = "climatology_"+season+".nc"
    setfileoption("nc","Format","LargeFile")
    if (fileexists(out_file)) then
    system("rm -f "+out_file)
    end if
    f_out = addfile(out_file,"c")
    f_out->track_density = track_density

