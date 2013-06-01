#! /bin/bash

for((i=0;i<32;i=i+1))
do
    echo -e "1\n$[$i+1]\n" | cat >proc_$i.inp
done
