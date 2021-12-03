#!/bin/sh  
exe="InfoScan" #你需要发布的程序名称
des="/public/home/meisq/05QT/08software/InfoScan" #创建文件夹的位置
deplist=$(ldd $exe | awk  '{if (match($3,"/")){ printf("%s "),$3 } }')  
cp $deplist $des
