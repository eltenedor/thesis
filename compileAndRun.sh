#! /bin/bash

currDir=$PWD

if [ -z "$1" ]
then
    iSt=2
    iEn=256
else
    iSt=$1
    if [ -z "$2" ]
    then
        iEn=$iSt
    else
        iEn=$2
    fi
fi


expath=../order_verification/test_convective

if [ ! -d $expath ]
then
    echo -e "creating $expath"
    mkdir $expath
fi

rm -f $expath/ERR.out

for ((i=$iSt;i<=$iEn;i=i*2))
do

    if [ ! -d $expath/$i ]
    then
        echo -e "creating $expath/$i"
        mkdir $expath/$i
    fi

    rm -f $expath/$i/grid.input
    echo -e "compiling and linking $i/.grgen"
    make grgen expath=$expath/$i >/dev/null
    cd $expath/$i
    echo -e "*\ngrid.out\n-1 1 $i\n-1 1 $i\n1 1 1 1\n" | cat>grid.input
    echo -e "running $i/.grgen"
    ./grgen <grid.input >/dev/null
    cd $currDir
    
    echo -e "compiling and linking $i/.solver"
    make solver expath=$expath/$i >/dev/null 
    cd $expath/$i
    echo -e "running $i/.solver"
    ./solver 1>/dev/null 2>log 
    echo -e "$i*$i: $(cat ERR.out)" | cat >>../ERR.out
    cd $currDir
done

#gedit $expath/ERR.out &

