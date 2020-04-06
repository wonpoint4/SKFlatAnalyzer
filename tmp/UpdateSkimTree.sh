#!/bin/bash

if [ -z "$1" ]
then 
    echo "Usage: $0 SKIMDIR [dry]"
    echo "Example: $0 /gv0/DATA/SKFlat/Run2Legacy_v3/2016/DATA_SkimTree_SMP"
    exit 1
fi

SKIMDIR=$(realpath $1)
if [ ! -z "$2" ]
then
    echo "Dry-run..."
    DRY=true
fi

while read line
do
    linesplit=(${line//\// })
    for i in $(seq 0 "${#linesplit[@]}")
    do
	if [ "${linesplit[$i]}" = SKFlat ]
	then
	    array=(${linesplit[@]:$i})
	    break
	fi
    done
    VERSION=${array[1]}
    YEAR=${array[2]}
    TYPE=$(echo ${array[3]}|grep -o "^[A-Z]*")
    SKIM=${array[3]#${TYPE}_}
    SAMPLE=${array[4]}
    if [ "$TYPE" = DATA ]
    then 
	[ ${#array[@]} -eq 7 ] || continue
	PERIOD=${array[5]}
	TAG=${array[6]}
	TARGET=$SKFlat_WD/data/$VERSION/$YEAR/Sample/ForSNU/${SKIM}_${SAMPLE}_${PERIOD#period}.txt
    elif [ "$TYPE" = MC ]
    then 
	[ ${#array[@]} -eq 6 ] || continue
	TAG=${array[5]}
	SAMPLE=$(grep $SAMPLE $SKFlat_WD/data/$VERSION/$YEAR/Sample/SampleSummary_MC.txt|awk '{print $1}'|head -n1)
	TARGET=$SKFlat_WD/data/$VERSION/$YEAR/Sample/ForSNU/${SKIM}_${SAMPLE}.txt
	[ "$SAMPLE" = "" ] && { echo "Cannot find alias for ${array[4]}"; continue; }
    fi	    
#    echo VERSION=$VERSION
#    echo YEAR=$YEAR
#    echo TYPE=$TYPE
#    echo SKIM=$SKIM
#    echo SAMPLE=$SAMPLE
#    echo PERIOD=$PERIOD
#    echo TAG=$TAG
    echo "find $line -type f |sort -V > $TARGET"
    [ "$DRY" != true ] && find $line -type f |sort -V > $TARGET
done< <(find $SKIMDIR -type d | sort -V )
