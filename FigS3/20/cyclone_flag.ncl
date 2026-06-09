;根据给定区域内的气旋格点数确定该区域是否受气旋影响
;若某日该区域内气旋格点数大于总格点数百分比阈值，则认为该区域受气旋影响，flag置1，否则置0
;输出文件cyclone_flag.nc，维度为(year,day)，year为1970-2022

year_start = 1970
year_end = 2022
total_year = ispan(year_start, year_end, 1)
season = "SON"

lon1 = 150
lon2 = 170
lat1 = -45
lat2 = -30
date_start = 817
date_end = 1215

ratio = 0.2 ;区域格点数阈值比例
degree = 0.5
threshold = ceil(((lon2-lon1+1)*(1/degree)*(lat2-lat1+1)*(1/degree))*ratio)

nyear = dimsizes(total_year)
nday = 121
flag = new((/nyear,nday/),"integer")
do year = year_start, year_end
    f = addfile("/g/data/if69/ls3248/data/weather/cyclone/"+season+"/cyclone."+year+"."+season+".nc","r")
    flag_temp = f->flag({year*10000+date_start:year*10000+date_end},{lat1:lat2},{lon1:lon2})

    ;逐日判断区域内气旋格点数
    do i = 0, nday-1
        flag_sum = sum(flag_temp(i,:,:))
        if (flag_sum .gt. threshold) then
            flag(year-year_start,i) = 1
        else
            flag(year-year_start,i) = 0
        end if
    end do
    print(year+" done")
end do

flag!0 = "year"
flag!1 = "day"
flag&year = total_year
flag&day = flag_temp&time - year_end*10000

;-----output the Comp_avg_hgt
out_file = "cyclone_flag.nc"

if (fileexists(out_file)) then
    system("rm -f "+out_file)
end if
f_out = addfile(out_file,"c")
f_out->flag = flag
