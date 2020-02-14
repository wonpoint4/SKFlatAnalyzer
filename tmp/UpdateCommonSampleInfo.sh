#!/bin/bash
while read line <&3
do
    array=(${line//\// })
    YEAR=${array[6]}
    NAME=$(echo ${array[7]/.root/}|sed 's/GetEffLumi_//')
    DASNAME=$(head -n1 $SKFlat_WD/data/$SKFlatV/$YEAR/Sample/ForSNU/${NAME}.txt|sed 's@/@ @g'|awk '{print $7}')
    CROSSSECTIONS=($(find $SKFlat_WD -name SampleSummary*.txt|xargs -i cat {}|grep "$NAME "|awk '{print $3}'|grep -v FIXMECROSSSECTION|uniq))
    echo YEAR=$YEAR NAME=$NAME DASNAME=$DASNAME CROSSSECTIONS=${CROSSSECTIONS[@]}
    if [ ${#CROSSSECTIONS[@]} -eq 1 ];
    then 
	CROSSSECTION=${CROSSSECTIONS[0]}
    else
	echo "candidate cross sections= ${CROSSSECTIONS[@]}"
	read -p "select cross section: " CROSSSECTION
	[ -z "$CROSSSECTION" ] && CROSSSECTION=FIXMECROSSSECTION
    fi
    NUMS=($(echo 'cout<<Form("%f\t%f\t",sumW->GetEntries(),sumW->GetSum())<<endl;'|root -b -l ${line}|tail -n1))
    OUT=$SKFlat_WD/data/$SKFlatV/$YEAR/Sample/CommonSampleInfo/${NAME}.txt
    echo " > Update $OUT with $line"
    echo $NAME $DASNAME $CROSSSECTION ${NUMS[0]} ${NUMS[1]}
    read -p "(y/n): " YES
    [ "$YES" = "y" ] || [ "$YES" = "Y" ] || { echo "terminated by user"; exit 1; }
    echo "# alias PD xsec nmc sumw" > $OUT
    echo $NAME $DASNAME $CROSSSECTION ${NUMS[0]} ${NUMS[1]} >> $OUT
    SUMMARYFILE=$SKFlat_WD/data/$SKFlatV/$YEAR/Sample/SampleSummary_MC.txt
    EXISTSUMMARY=$(grep $DASNAME $SUMMARYFILE|grep $NAME)
    if [[ $EXISTSUMMARY == *"FIXMECROSSSECTION"* ]]
    then
	echo "exist at summary with FIXMECROSSSECTION,"
	echo "delete lines"
	sed -i "/$DASSNAME.*FIXMECROSSSECTION/d" $SUMMARYFILE
	echo "add lines"
	echo $NAME $DASNAME $CROSSSECTION ${NUMS[0]} ${NUMS[1]}>>$SUMMARYFILE
    elif [ -z "$EXISTSUMMARY" ]
    then
	echo "no lines for summary"
	echo "add lines"
	echo $NAME $DASNAME $CROSSSECTION ${NUMS[0]} ${NUMS[1]}>>$SUMMARYFILE
    fi
done 3< <(find $SKFlatOutputDir$SKFlatV/GetEffLumi -type f|sort)



#find data/Run2Legacy_v3/*/Sample/CommonSampleInfo -type f -mtime -1|while read line;do array=(${line//\// });NAME=${array[5]/.txt/}; DASNAME=$(head -n1 ${line/CommonSampleInfo/ForSNU/}|awk 'BEGIN{FS="/"}{print $9}'); NUMS=($(echo 'cout<<sumW->GetSum()<<"\t"<<sumW->GetEntries()<<"\t"'|root -b -l $SKFlatOutputDir/$SKFlatV/GetEffLumi/${array[2]}/GetEffLumi_${array[5]/.txt/.root}|tail -n1|awk '{print $1"\t"$2}')); echo $NAME $DASNAME CROSSSECTION ${NUMS[@]} > $line;done
