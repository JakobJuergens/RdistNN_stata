version 17
do "tests/test_framework.do"
clear
set obs 12
gen double fr_x=_n/4
gen double fr_y=sin(fr_x)+fr_x
gen double fr_den=1+fr_x/10
matrix FR_NEW=(.25\1.25\2.75)
