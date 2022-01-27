#!/bin/sh

#RmATom='B'
#RmATomPosi="0.500000000    0.547626139    0.799401197"

#LISTRange='5 12 20 28 39'
#LISTHight='1 2 3'

#Directory='Bubble_z_'

#POSCARFILE='POSCAR_bubble_'
#DeFectFILE='POSCAR_bubble_VB_'

#NOWPWD=`pwd`

#function of change atom number
function ChangATOM {
    
    DefectfileName=${1}
    RemoveAtom=${2}

    ntypDATA=`sed -n '6p' ${DefectfileName}`
    natDATA=`sed -n '7p' ${DefectfileName}`
        
    ntype=$(echo $ntypDATA | wc -w)
        
    NEWntypDATA=""
    NEWnatDATA=""
    for i in `seq 1 $ntype`
    do
        ONEAtom=$(echo $ntypDATA  | cut -f $i -d ' ')

        ONEnat=$(echo $natDATA | cut -f $i -d ' ')
        if [ "$ONEAtom" == "${RemoveAtom}" ]; then
            ONEnat=$((ONEnat-1))
        fi
        echo $ONEnat
        if [ $ONEnat -ne 0 ]; then
            NEWnatDATA+="${ONEnat} "
            NEWntypDATA+="${ONEAtom} "
        fi
    done
    sed -i "6s/.*/$NEWntypDATA/g" ${DefectfileName}
    sed -i "7s/.*/$NEWnatDATA/g" ${DefectfileName}
}

#function for checing overlap atom
function CheckingList {
    List1="$1"
    List2="$2"
    overLap=""
    for word in $List1
    do
        if [[ "$List2" == *"$word"* ]];then
            overLap+="$word "
        fi
    done
    echo $overLap
}

#function of remove atom data only one
function FINDTargetATOMDATA {

    DefectfileName=${1}
    RmATOMX=${2}
    RmATOMY=${3}
    RmATOMZ=${4}
    
    ListX=$(cat $DefectfileName | awk "\$1 ~ /${RmATOMX}/ {print NR}")
    ListY=$(cat $DefectfileName | awk "\$2 ~ /${RmATOMY}/ {print NR}")
    ListZ=$(cat $DefectfileName | awk "\$3 ~ /${RmATOMZ}/ {print NR}")

    CandidateX=$(CheckingList "${ListX}" "$ListY")
    CandidateY=$(CheckingList "$ListY" "$ListX")
    Candidate=$(CheckingList "$CandidateX" "$CandidateY")

    NumCandi=$(echo $Candidate | wc -w)
    if (( $NumCandi == 0 ));then
        echo "There is no matched atom in POSCAR-FILE"
        echo "plz, change remove target atom position!!"
        echo "(X,Y,Z) : $RmATOMX, $RmATOMY, $RmATOMZ"
        exit 1
    elif (( $NumCandi == 1 ));then
        echo $Candidate
    else
        for word in $Candidate
        do
            if [[ "$ListZ" == *"$word"* ]];then
                echo $word
            fi
        done
    fi
}

function RemoveATOMDATA {
    DefectfileName=${1}
    Removeline=${2}

    RMATOMDATA=`sed -n ${Removeline}p $DefectfileName`

    echo "we remove ${Removeline}-line"
    echo "(X,Y,Z) :$RMATOMDATA"
    sed -i "${Removeline}d" $DefectfileName
}

function RMnewLINE {
    DefectfileName=${1}
    sed -i "s///g" $DefectfileName
}

