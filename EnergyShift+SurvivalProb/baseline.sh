#! /usr/bin/env bash

RUNNUMBER=0
for((i=0; i<2601; i++))
do


        ds=`echo $((  ${RUNNUMBER}  /  51  ))` #변수 정의 할 때 = 사이에 띄어 쓰기 있으면 안 됨
        dm=`echo $((  ${RUNNUMBER}  %  51  ))`

        RUNNUMBER=`echo $((RUNNUMBER+1))`

        root -l -q getBinContent_baseline2.C `echo $(( ${ds} ))` `echo $(( ${dm} ))`
        echo `echo $(( ${ds} ))` `echo $(( ${dm} ))`
       
done

