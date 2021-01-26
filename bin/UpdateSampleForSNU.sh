#!/bin/bash
if [ -z "$SKFlatV" ];then
    echo "Please source setup.sh"
    exit 1
fi

if [ -z "$1" ]; then
    echo "Usage: $0 SEARCHDIR"
    echo "Example: $0 /gv0/DATA/SKFlat/$SKFlatV"
    exit 1
fi

SEARCHDIR=$(realpath $1)
if [ ! -z "$2" ];then
    echo "Dry-run..."
    DRY=true
fi

## find /gv0/DATA/SKFlat/Run2UltraLegacy_v1/2017/DATA/SingleElectron/periodB/210120_042549/ -type f |sort -V > $SKFlat_WD/data/Run2UltraLegacy_v1/2017/Sample/ForSNU/
while read line <&3; do
    array=(${line//\// })
    for i in $(seq "${#array[@]}"); do
	if [ "${array[$i]}" = "$SKFlatV" ]; then
	    VI=$i ## SKFlatV position index
	    break
	fi
    done
    YEAR="${array[$(($VI+1))]}"
    TYPE="${array[$(($VI+2))]}"
    if [[ "$TYPE" = *"DATA"* ]]; then
	test "${#array[@]}" -eq "$(($VI+5))" || continue
	SAMPLE="${array[$(($VI+3))]}"
	PERIOD="${array[$(($VI+4))]}"
	ALIAS="${SAMPLE}_${PERIOD//period/}"
    elif [[ "$TYPE" = *"MC"* ]]; then
	test "${#array[@]}" -eq "$(($VI+4))" || continue
	SAMPLE="${array[$(($VI+3))]}"
	ALIAS=$(grep $SAMPLE $SKFlat_WD/data/$SKFlatV/$YEAR/Sample/SampleSummary_MC.txt|awk '{print $1}'|head -n1)
    else continue;
    fi
    CRABTAG=$(ls $line|sort -V|tail -n1)
    test -z $CRABTAG && continue
    SOURCE=$line/$CRABTAG
    TARGET=$SKFlat_WD/data/$SKFlatV/$YEAR/Sample/ForSNU/$ALIAS.txt
    if [ -f $TARGET ];then
	diffout="$(diff <(find $SOURCE -type f|sort -V) $TARGET)"
	if [ $(echo "$diffout"|wc -l) -gt 1 ]; then
	    echo "$diffout"|tail
	else continue
	fi
    fi
    echo "find $SOURCE -type f|sort -V > $SKFlat_WD/data/$SKFlatV/$YEAR/Sample/ForSNU/$ALIAS.txt"
    if [ "$DRY" != "true" ]; then
	read -p "exec? (y/n): " YES
	if [ "$YES" = "y" ];then
	    find $SOURCE -type f|sort -V > $SKFlat_WD/data/$SKFlatV/$YEAR/Sample/ForSNU/$ALIAS.txt
	fi
    fi
done 3< <(find $SEARCHDIR -type d|grep $SKFlatV|sort -V)
