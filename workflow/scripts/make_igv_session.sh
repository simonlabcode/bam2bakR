#!/bin/bash
source $1

if [[ -n $GNUPARALLEL_MOD ]]; then module load ${GNUPARALLEL_MOD}; fi

case $SPECIES in
    Hs) genome="hg38" ;;
    Mm) genome="mm10" ;;
    Dm) genome="dm6" ;;
    *)  genome="hg38" ;;
esac

muts=$(echo $mut_tracks | tr ',' ' ')

cd $MASTER_DIR/tracks.dir

echo '<?xml version="1.0" encoding="UTF-8" standalone="no"?>' > igv_session.xml
echo '<Session genome="'$genome'" hasGeneTrack="true" hasSequenceTrack="true" locus="All" nextAutoscaleGroup="3" version="8">' >> igv_session.xml



echo '    <Resources>' >> igv_session.xml

for sample in $(ls *.tdf); do
    echo '        <Resource path="'$sample'"/>' >> igv_session.xml
done

echo '    </Resources>' >> igv_session.xml



echo '    <Panel height="482" name="DataPanel" width="1363">' >> igv_session.xml

for sample in ${samples[@]}; do

    for b in $muts; do
        if [ $b == "GA" ]; then
            colVal[0]='200,200,200'
            colVal[1]='120,188,230'
            colVal[2]='65,125,195'
            colVal[3]='36,110,182'
            colVal[4]='27,78,165'
            colVal[5]='18,50,120'
        else
            colVal[0]='200,200,200'
            colVal[1]='250,150,150'
            colVal[2]='250,0,0'
            colVal[3]='150,0,0'
            colVal[4]='100,0,0'
            colVal[5]='50,0,0'
        fi

        echo '        <Track attributeKey="Overlay" autoScale="false" autoscaleGroup="1" clazz="org.broad.igv.track.MergedTracks" fontSize="10" height="80" id="id'$sample'" name="'$sample' '$b'" renderer="BAR_CHART" visible="true">' >> igv_session.xml

        $GNUPARALLEL -j 1 "echo '            <Track {5}=\"{6}\" attributeKey=\"{1} {2} {7} {3}\" autoScale=\"false\" autoscaleGroup=\"1\" fontSize=\"10\" id=\"{1}.{2}.{7}.{4}.tdf\" name=\"{1} {2} {7} {3}\" renderer=\"BAR_CHART\" visible=\"true\" windowFunction=\"mean\"/>'" \
                    ::: $sample \
                    ::: $b \
                    ::: plus minus \
                    :::+ pos min \
                    :::+ color altColor \
                    ::: ${colVal[@]} \
                    :::+ $(seq 0 5) >> igv_session.xml

        echo '        </Track>' >> igv_session.xml
    done

done

# For Mutation position tracks
if [ $mut_pos = "TRUE" ]; then
    for sample in ${samples[@]}; do

        for b in $muts; do
            if [ $b == "GA" ]; then
                colVal='65,125,195'
            else
                colVal='250,0,0'
            fi

            echo '        <Track attributeKey="Overlay" autoScale="false" autoscaleGroup="2" clazz="org.broad.igv.track.MergedTracks" fontSize="10" height="40" id="idmut'$sample'" name="Mut position '$sample' '$b'" renderer="BAR_CHART" visible="true">' >> igv_session.xml
            echo '            <Track color="'${colVal}'" attributeKey="mutPos '${sample}' '${b}' plus" autoScale="false" autoscaleGroup="2" fontSize="10" id="'${sample}'.'${b}'.muts.pos.tdf" name="Mut position '${sample}' '${b}' plus" renderer="BAR_CHART" visible="true" windowFunction="mean"/>' >> igv_session.xml
            echo '            <Track altColor="'${colVal}'" attributeKey="mutPos '${sample}' '${b}' minus" autoScale="false" autoscaleGroup="2" fontSize="10" id="'${sample}'.'${b}'.muts.min.tdf" name="Mut position '${sample}' '${b}' minus" renderer="BAR_CHART" visible="true" windowFunction="mean"/>' >> igv_session.xml
            echo '        </Track>' >> igv_session.xml
        done

    done
fi

echo '    </Panel>' >> igv_session.xml



echo '    <Panel height="238" name="FeaturePanel" width="1363">' >> igv_session.xml
echo '        <Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="10" id="Reference sequence" name="Reference sequence" visible="true"/>' >> igv_session.xml
echo '        <Track altColor="0,0,178" attributeKey="Gene" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;845.0;255,255,255;0,0,178" fontSize="10" height="35" id="'$genome'_genes" name="Gene" visible="true"/>' >> igv_session.xml
echo '    </Panel>' >> igv_session.xml

echo '    <PanelLayout dividerFractions="0.6657496561210454"/>' >> igv_session.xml
echo '</Session>' >> igv_session.xml

