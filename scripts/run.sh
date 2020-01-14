#!/bin/bash
BM_PATH=${BENCHMARK_PATH?"Need to set environment variable BENCHMARK_PATH"}

BM_ABBR=(f{01..12})
BM_LONG=(FPGA{01..12})
# BM_ABBR=(e{1..4} f{01..12} c{1..5} cf{01..13})
# BM_LONG=(FPGA-example{1..4} FPGA{01..12} clk_design{1..5} CLK-FPGA{01..13})
BM16_NUM=16
BM_NUM=${#BM_ABBR[@]}

if [ -z $1 ] ; then
    echo "Usage:"
    echo "$0 [mode] <benchmark>|all [option1] [option2]"
    echo "Mode:"
    echo " tdm_part     partition the circuit and gen sub problems"
    echo " tdm_place    solve sub problems by global placement"
    echo " tdm_time     solve tdm optimization"
    echo " gdb          debug with gdb"
    echo " valgrind     check memory error with valgrind"
    echo " vgdb         run valgrind with gdb debugging"
    echo "Options:"
    echo "-flow         select flow"
    echo "-lg           select discretization method"
    echo "-thread       select number of threads"
    echo "Available benchmarks:"
    for ((i=0; i<$BM_NUM; ++i)); do
        echo ${BM_ABBR[i]} - ${BM_LONG[i]}
    done
    exit
fi

PREFIX=""
if [ $1 = "gdb" ] ; then
    shift
    PREFIX="gdb --args "
    MODE=$1_gdb
    shift
elif [ $1 = "valgrind" ] ; then
    shift
    PREFIX="valgrind "
    MODE=$1_valgrind
    shift
elif [ $1 = "vgdb" ] ; then
    shift
    PREFIX="valgrind --vgdb-error=0"
    MODE=$1_vgdb
    shift
elif [ $1 = "tdm_place" -o $1 = "tdm_time" -o $1 = "tdm_part" ] ; then
    MODE=$1
    shift
else
    echo "unknown mode"
    exit
fi
echo "Mode: $MODE"

DATE=`date +"%m%d"`
BM=$1
shift
if [[ $MODE == *"tdm_place"* || $MODE == *"tdm_time"* || $MODE == *"tdm_part"* ]] ; then 
    N_SUBPROBLEM=$1
    if [ -z $1 ] ; then
        echo "need to enter the number of devices"
        exit
    fi
    shift
fi
OPTIONS="$@"
BINARY="larf"

FOUND="false"

for ((i=0; i<$BM_NUM; ++i)); do
    if [ $BM = "all" -o $BM = ${BM_ABBR[i]} \
    -o \( $BM = "all16" -a $i -lt $BM16_NUM \) \
    -o \( $BM = "all17" -a $i -ge $BM16_NUM \) ] ; then
        if [ $i -lt $BM16_NUM ] ; then
            YEAR=16
        else
            YEAR=17
        fi

        if [[ $MODE == *"tdm_place"* ]] ; then
            for ((SUB_PROBLEM=0; SUB_PROBLEM<$N_SUBPROBLEM; ++SUB_PROBLEM)); do
                BENCHMARK=${BM_ABBR[i]}     
                INPUT_AUX=${SUB_PROBLEM}/design.aux
                OUTPUT_PL=${BM_LONG[i]}_${SUB_PROBLEM}.pl    
                BM_DATED=${BENCHMARK}_${SUB_PROBLEM}_${DATE}

                cd $BENCHMARK/
                echo "$PREFIX ../$BINARY -aux $INPUT_AUX -out $OUTPUT_PL -flow tdm_place"
                    $PREFIX ../$BINARY -aux $INPUT_AUX -out $OUTPUT_PL -flow tdm_place | tee $BM_DATED.log
                cd ../
            done
        elif [[ $MODE == *"tdm_time"* ]] ; then
            BENCHMARK=${BM_ABBR[i]}
            cd $BENCHMARK/
            mkdir -p combine
            cp design_inter.nets combine/design.nets
            cp $BM_PATH/ispd20$YEAR/${BM_LONG[i]}/design.aux combine/design.aux
            cp $BM_PATH/ispd20$YEAR/${BM_LONG[i]}/design.nodes combine/design.nodes
            cp $BM_PATH/ispd20$YEAR/${BM_LONG[i]}/design.lib combine/design.lib
            cp $BM_PATH/ispd20$YEAR/${BM_LONG[i]}/design.wts combine/design.wts
            cp $BM_PATH/ispd20$YEAR/${BM_LONG[i]}/design.scl combine/design.scl
            cp $BM_PATH/ispd20$YEAR/${BM_LONG[i]}/design.pl combine/design.pl
            > instance.pos
            for ((SUB_PROBLEM=0; SUB_PROBLEM<$N_SUBPROBLEM; ++SUB_PROBLEM)); do
                cat $SUB_PROBLEM/design.nets >> combine/design.nets
                cat ${BM_LONG[i]}_${SUB_PROBLEM}.pl >> instance.pos
            done

            OUTPUT_FILE=${BM_ABBR[i]}.tdm
            INPUT_AUX=combine/design.aux
            BM_DATED=${BENCHMARK}_combine_${DATE}
            echo "$PREFIX ../$BINARY -aux $INPUT_AUX -flow tdm_time -out $OUTPUT_FILE -partition $N_SUBPROBLEM $OPTIONS"
                    $PREFIX ../$BINARY -aux $INPUT_AUX -flow tdm_time -out $OUTPUT_FILE -partition $N_SUBPROBLEM $OPTIONS | tee $BM_DATED.log
            cd ../
        elif [[ $MODE == *"tdm_part"* ]] ; then
            INPUT_AUX=$BM_PATH/ispd20$YEAR/${BM_LONG[i]}/design.aux
            BENCHMARK=${BM_ABBR[i]}
            BM_DATED=${BENCHMARK}_${DATE}

            mkdir -p $BENCHMARK/
            cd $BENCHMARK/
            echo "$PREFIX ../$BINARY -aux $INPUT_AUX -flow tdm_part -partition $N_SUBPROBLEM $OPTIONS"
                $PREFIX ../$BINARY -aux $INPUT_AUX -flow tdm_part -partition $N_SUBPROBLEM $OPTIONS| tee $BM_DATED.log
            cd ../
        fi
        
        FOUND="yes"
    fi
done

if [ $FOUND != "yes" ] ; then
    echo "ERROR: no benchmark with name: $BM"
fi
