;-----SON define
; EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016, 2024/)
; CP_ElNino_year = (/1987, 1992, 1995, 2003, 2005, 2010, 2019/)
; LaNina_year = (/1984, 1989, 1999, 2001, 2011, 2012, 2022, 2023/) ;CP
EP_ElNino_year = (/1983, 1988, 1998, 2007, 2016/)
CP_ElNino_year = (/1987, 1992, 1995, 2003, 2005, 2010, 2019/)
EP_LaNina_year = (/1996, 2000, 2008, 2021/)
CP_LaNina_year = (/1984, 1989, 1999, 2011, 2012, 2022, 2023/)
LaNina_year = array_append_record(EP_LaNina_year,CP_LaNina_year,0)

neutral_year = (/1981, 1982, 1986, 1990, 1991, 1993, 1994, 1997, 2002, 2004, 2013, 2014, 2020/)
nEPElNino = dimsizes(EP_ElNino_year)
nCPElNino = dimsizes(CP_ElNino_year)    
nLaNina = dimsizes(LaNina_year)
nneutral = dimsizes(neutral_year)

;-----standardrize the entire series SON-SON-...-SON
yearStart = 1980
yearEnd = 2024
years = ispan(yearStart,yearEnd,1)
day = ispan(1,121,1) ;17 Aug - SON - 15 Dec
nyear = dimsizes(years)
nday = dimsizes(day)

f = addfile("/g/data/m35/ls3248/ENSO_cyclone/paper/Fig2/pattern_index.nc","r") ;pattern correlation with 3 cluster patterns within all years from 1980-2024
;-----standardization  
EPElNino = f->p_cor_EPElNino
EPElNino_index = onedtond(dim_standardize(ndtooned(EPElNino), 0), (/nyear,nday/))
CPElNino = f->p_cor_CPElNino
CPElNino_index = onedtond(dim_standardize(ndtooned(CPElNino), 0), (/nyear,nday/))
LaNina = f->p_cor_LaNina
LaNina_index = onedtond(dim_standardize(ndtooned(LaNina), 0), (/nyear,nday/))
copy_VarCoords(EPElNino, EPElNino_index)
copy_VarCoords(CPElNino, CPElNino_index)
copy_VarCoords(LaNina, LaNina_index)

min_duration = 10  ; Minimum consecutive days
threshold = 1. ; threshold for pattern index

;----- for EPElNino_in_EPElNino_year
count = 0
count_all = new((/nEPElNino/),"integer")
count_all = 0
do i = 0,nEPElNino-1
    year = EP_ElNino_year(i)-1
    index = where(EPElNino_index({year},:).gt.threshold, 1, 0)
    nday = dimsizes(index)
    do j = 0,nday-1
        if index(j) .eq. 1 then
            count = count + 1
        end if
        if index(j) .eq. 0 then
            if count .ge. min_duration then
                count_all(i) = count_all(i) + count
            end if
            count = 0
        end if
    end do
    count = 0
end do
percentage_EPElNino_in_EPElNino_year = count_all/121.*100.
EPElNino_inEPElNino = avg(percentage_EPElNino_in_EPElNino_year)
EPElNino_inEPElNino_err = stddev(percentage_EPElNino_in_EPElNino_year)/sqrt(nEPElNino)
delete(count_all)

;----- for CPElNino_in_EPElNino_year
count = 0
count_all = new((/nEPElNino/),"integer")
count_all = 0
do i = 0,nEPElNino-1
    year = EP_ElNino_year(i)-1
    index = where(CPElNino_index({year},:).gt.threshold, 1, 0)
    nday = dimsizes(index)
    do j = 0,nday-1
        if index(j) .eq. 1 then
            count = count + 1
        end if
        if index(j) .eq. 0 then
            if count .ge. min_duration then
                count_all(i) = count_all(i) + count
            end if
            count = 0
        end if
    end do
    count = 0
end do
percentage_CPElNino_in_EPElNino_year = count_all/121.*100.
CPElNino_inEPElNino = avg(percentage_CPElNino_in_EPElNino_year)
CPElNino_inEPElNino_err = stddev(percentage_CPElNino_in_EPElNino_year)/sqrt(nEPElNino)
delete(count_all)

;----- for LaNina_in_EPElNino_year
count = 0
count_all = new((/nEPElNino/),"integer")
count_all = 0
do i = 0,nEPElNino-1
    year = EP_ElNino_year(i)-1
    index = where(LaNina_index({year},:).gt.threshold, 1, 0)
    nday = dimsizes(index)
    do j = 0,nday-1
        if index(j) .eq. 1 then
            count = count + 1
        end if
        if index(j) .eq. 0 then
            if count .ge. min_duration then
                count_all(i) = count_all(i) + count
            end if
            count = 0
        end if
    end do
    count = 0
end do
percentage_LaNina_in_EPElNino_year = count_all/121.*100.
LaNina_inEPElNino = avg(percentage_LaNina_in_EPElNino_year)
LaNina_inEPElNino_err = stddev(percentage_LaNina_in_EPElNino_year)/sqrt(nEPElNino)
delete(count_all)

;----- for EPElNino_in_CPElNino_year
count = 0
count_all = new((/nCPElNino/),"integer")
count_all = 0
do i = 0,nCPElNino-1
    year = CP_ElNino_year(i)-1
    index = where(EPElNino_index({year},:).gt.threshold, 1, 0)
    nday = dimsizes(index)
    do j = 0,nday-1
        if index(j) .eq. 1 then
            count = count + 1
        end if
        if index(j) .eq. 0 then
            if count .ge. min_duration then
                count_all(i) = count_all(i) + count
            end if
            count = 0
        end if
    end do
    count = 0
end do
percentage_EPElNino_in_CPElNino_year = count_all/121.*100.
EPElNino_inCPElNino = avg(percentage_EPElNino_in_CPElNino_year)
EPElNino_inCPElNino_err = stddev(percentage_EPElNino_in_CPElNino_year)/sqrt(nCPElNino)
delete(count_all)

;----- for CPElNino_in_CPElNino_year
count = 0
count_all = new((/nCPElNino/),"integer")
count_all = 0
do i = 0,nCPElNino-1
    year = CP_ElNino_year(i)-1
    index = where(CPElNino_index({year},:).gt.threshold, 1, 0)
    nday = dimsizes(index)
    do j = 0,nday-1
        if index(j) .eq. 1 then
            count = count + 1
        end if
        if index(j) .eq. 0 then
            if count .ge. min_duration then
                count_all(i) = count_all(i) + count
            end if
            count = 0
        end if
    end do
    count = 0
end do
percentage_CPElNino_in_CPElNino_year = count_all/121.*100.
CPElNino_inCPElNino = avg(percentage_CPElNino_in_CPElNino_year)
CPElNino_inCPElNino_err = stddev(percentage_CPElNino_in_CPElNino_year)/sqrt(nCPElNino)
delete(count_all)

;----- for LaNina_in_CPElNino_year
count = 0
count_all = new((/nCPElNino/),"integer")
count_all = 0
do i = 0,nCPElNino-1
    year = CP_ElNino_year(i)-1
    index = where(LaNina_index({year},:).gt.threshold, 1, 0)
    nday = dimsizes(index)
    do j = 0,nday-1
        if index(j) .eq. 1 then
            count = count + 1
        end if
        if index(j) .eq. 0 then
            if count .ge. min_duration then
                count_all(i) = count_all(i) + count
            end if
            count = 0
        end if
    end do
    count = 0
end do
percentage_LaNina_in_CPElNino_year = count_all/121.*100.
LaNina_inCPElNino = avg(percentage_LaNina_in_CPElNino_year)
LaNina_inCPElNino_err = stddev(percentage_LaNina_in_CPElNino_year)/sqrt(nCPElNino)
delete(count_all)

;----- for EPElNino_in_LaNina_year
count = 0
count_all = new((/nLaNina/),"integer")
count_all = 0
do i = 0,nLaNina-1
    year = LaNina_year(i)-1
    index = where(EPElNino_index({year},:).gt.threshold, 1, 0)
    nday = dimsizes(index)
    do j = 0,nday-1
        if index(j) .eq. 1 then
            count = count + 1
        end if
        if index(j) .eq. 0 then
            if count .ge. min_duration then
                count_all(i) = count_all(i) + count
            end if
            count = 0
        end if
    end do
    count = 0
end do
percentage_EPElNino_in_LaNina_year = count_all/121.*100.
EPElNino_inLaNina = avg(percentage_EPElNino_in_LaNina_year)
EPElNino_inLaNina_err = stddev(percentage_EPElNino_in_LaNina_year)/sqrt(nLaNina)
delete(count_all)

;----- for CPElNino_in_LaNina_year
count = 0
count_all = new((/nLaNina/),"integer")
count_all = 0
do i = 0,nLaNina-1
    year = LaNina_year(i)-1
    index = where(CPElNino_index({year},:).gt.threshold, 1, 0)
    nday = dimsizes(index)
    do j = 0,nday-1
        if index(j) .eq. 1 then
            count = count + 1
        end if
        if index(j) .eq. 0 then
            if count .ge. min_duration then
                count_all(i) = count_all(i) + count
            end if
            count = 0
        end if
    end do
    count = 0
end do
percentage_CPElNino_in_LaNina_year = count_all/121.*100.
CPElNino_inLaNina = avg(percentage_CPElNino_in_LaNina_year)
CPElNino_inLaNina_err = stddev(percentage_CPElNino_in_LaNina_year)/sqrt(nLaNina)
delete(count_all)

;----- for LaNina_in_LaNina_year
count = 0
count_all = new((/nLaNina/),"integer")
count_all = 0
do i = 0,nLaNina-1
    year = LaNina_year(i)-1
    index = where(LaNina_index({year},:).gt.threshold, 1, 0)
    nday = dimsizes(index)
    do j = 0,nday-1
        if index(j) .eq. 1 then
            count = count + 1
        end if
        if index(j) .eq. 0 then
            if count .ge. min_duration then
                count_all(i) = count_all(i) + count
            end if
            count = 0
        end if
    end do
    count = 0
end do
percentage_LaNina_in_LaNina_year = count_all/121.*100.
LaNina_inLaNina = avg(percentage_LaNina_in_LaNina_year)
LaNina_inLaNina_err = stddev(percentage_LaNina_in_LaNina_year)/sqrt(nLaNina)
delete(count_all)

;----- for EPElNino_in_Neutral_year
count = 0
count_all = new((/nneutral/),"integer")
count_all = 0
do i = 0,nneutral-1
    year = neutral_year(i)-1
    index = where(EPElNino_index({year},:).gt.threshold, 1, 0)
    nday = dimsizes(index)
    do j = 0,nday-1
        if index(j) .eq. 1 then
            count = count + 1
        end if
        if index(j) .eq. 0 then
            if count .ge. min_duration then
                count_all(i) = count_all(i) + count
            end if
            count = 0
        end if
    end do
    count = 0
end do
percentage_EPElNino_in_neutral_year = count_all/121.*100.
EPElNino_inNeutral = avg(percentage_EPElNino_in_neutral_year)
EPElNino_inNeutral_err = stddev(percentage_EPElNino_in_neutral_year)/sqrt(nneutral)
delete(count_all)

;----- for CPElNino_in_Neutral_year
count = 0
count_all = new((/nneutral/),"integer")
count_all = 0
do i = 0,nneutral-1
    year = neutral_year(i)-1
    index = where(CPElNino_index({year},:).gt.threshold, 1, 0)
    nday = dimsizes(index)
    do j = 0,nday-1
        if index(j) .eq. 1 then
            count = count + 1
        end if
        if index(j) .eq. 0 then
            if count .ge. min_duration then
                count_all(i) = count_all(i) + count
            end if
            count = 0
        end if
    end do
    count = 0
end do
percentage_CPElNino_in_neutral_year = count_all/121.*100.
CPElNino_inNeutral = avg(percentage_CPElNino_in_neutral_year)
CPElNino_inNeutral_err = stddev(percentage_CPElNino_in_neutral_year)/sqrt(nneutral)
delete(count_all)

;----- for LaNina_in_Neutral_year
count = 0
count_all = new((/nneutral/),"integer")
count_all = 0
do i = 0,nneutral-1
    year = neutral_year(i)-1
    index = where(LaNina_index({year},:).gt.threshold, 1, 0)
    nday = dimsizes(index)
    do j = 0,nday-1
        if index(j) .eq. 1 then
            count = count + 1
        end if
        if index(j) .eq. 0 then
            if count .ge. min_duration then
                count_all(i) = count_all(i) + count
            end if
            count = 0
        end if
    end do
    count = 0
end do
percentage_LaNina_in_neutral_year = count_all/121.*100.
LaNina_inNeutral = avg(percentage_LaNina_in_neutral_year)
LaNina_inNeutral_err = stddev(percentage_LaNina_in_neutral_year)/sqrt(nneutral)
delete(count_all)

print("EP pattern")
print("prob: "+EPElNino_inEPElNino+","+EPElNino_inCPElNino+","+EPElNino_inLaNina+","+EPElNino_inNeutral)

print("CP pattern")
print("prob: "+CPElNino_inEPElNino+","+CPElNino_inCPElNino+","+CPElNino_inLaNina+","+CPElNino_inNeutral)

print("LaNina pattern")
print("prob: "+LaNina_inEPElNino+","+LaNina_inCPElNino+","+LaNina_inLaNina+","+LaNina_inNeutral)

print("EP err")
print("prob: "+EPElNino_inEPElNino_err+","+EPElNino_inCPElNino_err+","+EPElNino_inLaNina_err+","+EPElNino_inNeutral_err)

print("CP err")
print("prob: "+CPElNino_inEPElNino_err+","+CPElNino_inCPElNino_err+","+CPElNino_inLaNina_err+","+CPElNino_inNeutral_err)

print("LaNina err")
print("prob: "+LaNina_inEPElNino_err+","+LaNina_inCPElNino_err+","+LaNina_inLaNina_err+","+LaNina_inNeutral_err)