#!/bin/sh

#RmATom='B'
#RmATomPosi="0.500000000    0.547626139    0.799401197"

#ChATom='N'
#limitRange='0.05' #Not angstrom, fraction
#NewATom='n'

#Directory='Bubble_z_'

#POSCARFILE='POSCAR_bubble_28'
#DeFectFILE='POSCAR_bubble_VB_28'


#####

#search New atom candidate
function SearchNewATOMDATA {
    DefectfileName=${1}
    OriginPosiX=${2}
    OriginPosiY=${3}
    OriginPosiZ=${4}

    limitRange=${5}

    
    LastLine=`cat $DefectfileName | wc -l`

    Candidate=""
    for i in `seq 9 $LastLine`
    do
        posi=`sed -n "$i"p < $DefectfileName`
        posiX=$(echo $(echo $posi| cut -f 1 -d ' ') - $OriginPosiX |bc)
        posiY=$(echo $(echo $posi| cut -f 2 -d ' ') - $OriginPosiY |bc)
        Distance=$(printf "%0.7f" $(echo "sqrt($posiX^2+$posiY^2)" | bc))
        
        if (( $(echo  "${Distance} < ${limitRange}" | bc -l) ));then
            Candidate+="${i} "
        fi

    done

    echo $Candidate
}


#change new atom data
function ChangeNewATOMDATA {
    DefectfileName=${1}
    ChangeAtom=${2}
    NewAtom=${3}
    LineNewAtom=${4}
    
    NumNewAtom=$(echo $LineNewAtom | wc -w)

    ntypDATA=`sed -n '6p' ${DefectfileName}`
    natDATA=`sed -n '7p' ${DefectfileName}`
        
    ntype=$(echo $ntypDATA | wc -w)
        
    #about nat, ntyp
    NEWntypDATA=""
    NEWnatDATA=""
    for i in `seq 1 $ntype`
    do
        ONEAtom=$(echo $ntypDATA  | cut -f $i -d ' ')

        ONEnat=$(echo $natDATA | cut -f $i -d ' ')
        if [ "$ONEAtom" == "${ChangeAtom}" ]; then
            ONEnat=$((ONEnat-$NumNewAtom))
        fi
        if [ $ONEnat -ne 0 ]; then
            NEWnatDATA+="${ONEnat} "
            NEWntypDATA+="${ONEAtom} "
        fi
    done
    NEWnatDATA+="${NumNewAtom} "
    NEWntypDATA+="${NewAtom} "

    sed -i "6s/.*/$NEWntypDATA/g" ${DefectfileName}
    sed -i "7s/.*/$NEWnatDATA/g" ${DefectfileName}

    #about atomic position
    reverseLine=$(echo $LineNewAtom | tr " " "\n" | sort -gr)
    for i in $reverseLine
    do
        Content=$(sed -n "$i"p $DefectfileName)
        echo $Content >> $DefectfileName
        sed -i "${i}d" $DefectfileName
    done
}

#####
