yearStart = 1980
yearEnd = 2024
totalYear = ispan(yearStart,yearEnd,1)

latN = -20
latS = -80
lonL = 80 ;80E
lonR = 280 ;80w

nyear = dimsizes(totalYear)
nday = 121
p_cor = new((/nyear,nday/),"float")
p_cor!0 = "year"
p_cor!1 = "day"
p_cor&year = totalYear
p_cor&day = ispan(1,121,1) ;from 8.17 to 12.15

opt = 1 ;1-uncentered, 0-centered, for anomaly pattern correlation, use 1

f = addfile("./composite_pattern.nc","r")
pattern = f->composite_pattern(:,::-1,:)

k = 0
do year = yearStart,yearEnd

    f = addfile("/g/data/if69/ls3248/data/ERA5/daily/z/LF/z.ERA5.SON.anom.LF."+year+".nc", "r")
    z = f->z(:,{500},:,:)
    hgt_temp = z/9.8
    copy_VarCoords(z, hgt_temp)

    ;-----interpolation
    LON = ispan(0,359,1)
    LAT = ispan(-90,90,1)
    LON@units = "degree_east"
    LAT@units = "degree_north"
    hgt = linint2_Wrap(hgt_temp&lon,hgt_temp&lat,hgt_temp, True, LON,LAT, 0)

    x = hgt(:,{latS:latN},{lonL:lonR})
    y = pattern(0,{latS:latN},{lonL:lonR}) ;0-EPElNino, 1-CPElNino, 2-LaNina


    lat = y&latitude
    pi = 3.1415926
    w = cos(pi/180.*lat)
    p_cor(k,:) = tofloat(pattern_cor(x,conform(x, y, (/1,2/)),w,opt))

    print(year+" is ok.")
    k = k+1
end do
p_cor_EPElNino = p_cor

k = 0
do year = yearStart,yearEnd

    f = addfile("/g/data/if69/ls3248/data/ERA5/daily/z/LF/z.ERA5.SON.anom.LF."+year+".nc", "r")
    z = f->z(:,{500},:,:)
    hgt_temp = z/9.8
    copy_VarCoords(z, hgt_temp)

    ;-----interpolation
    LON = ispan(0,359,1)
    LAT = ispan(-90,90,1)
    LON@units = "degree_east"
    LAT@units = "degree_north"
    hgt = linint2_Wrap(hgt_temp&lon,hgt_temp&lat,hgt_temp, True, LON,LAT, 0)

    x = hgt(:,{latS:latN},{lonL:lonR})
    y = pattern(1,{latS:latN},{lonL:lonR}) ;0-EPElNino, 1-CPElNino, 2-LaNina
    

    lat = y&latitude
    pi = 3.1415926
    w = cos(pi/180.*lat)
    p_cor(k,:) = tofloat(pattern_cor(x,conform(x, y, (/1,2/)),w,opt))

    print(year+" is ok.")
    k = k+1
end do
p_cor_CPElNino = p_cor

k = 0
do year = yearStart,yearEnd

    f = addfile("/g/data/if69/ls3248/data/ERA5/daily/z/LF/z.ERA5.SON.anom.LF."+year+".nc", "r")
    z = f->z(:,{500},:,:)
    hgt_temp = z/9.8
    copy_VarCoords(z, hgt_temp)

    ;-----interpolation
    LON = ispan(0,359,1)
    LAT = ispan(-90,90,1)
    LON@units = "degree_east"
    LAT@units = "degree_north"
    hgt = linint2_Wrap(hgt_temp&lon,hgt_temp&lat,hgt_temp, True, LON,LAT, 0)

    x = hgt(:,{latS:latN},{lonL:lonR})
    y = pattern(2,{latS:latN},{lonL:lonR}) ;0-EPElNino, 1-CPElNino, 2-LaNina
    

    lat = y&latitude
    pi = 3.1415926
    w = cos(pi/180.*lat)
    p_cor(k,:) = tofloat(pattern_cor(x,conform(x, y, (/1,2/)),w,opt))

    print(year+" is ok.")
    k = k+1
end do
p_cor_LaNina = p_cor

;-----output the Comp_avg_hgt
out_file = "pattern_index.nc"
if (fileexists(out_file)) then
    system("rm -f "+out_file)
end if
f_out = addfile(out_file,"c")
f_out->p_cor_EPElNino = p_cor_EPElNino
f_out->p_cor_CPElNino = p_cor_CPElNino
f_out->p_cor_LaNina = p_cor_LaNina

