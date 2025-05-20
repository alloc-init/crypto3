#!/bin/bash

timelimit="5m"

if [[ $# -ne 1 ]]; then
    echo "ERROR: csv output file required"
    exit 1
fi

csv=$1

if [[ $(uname) == Darwin ]]; then
    arch=mac
else
    arch=linux
fi

if ! which timeout 2>&1>/dev/null; then
    echo "ERROR: missing required command timeout (brew install coreutils)"
    exit 1
fi

if ! which fd 2>&1>/dev/null; then
    echo "ERROR: missing required command fd (brew install fd)"
    exit 1
fi

function list_make_targets() {
    make -qp \
        | awk -F':' '/^[a-zA-Z0-9][^$#\/\t=]*:([^=]|$)/ {split($1,A,/ /);for(i in A)print A[i]}' \
        | sort -u
}

raw_tests=$(list_make_targets | grep _test | sort)

# drop tests which we've logged already
if [[ $(wc -l "$csv" | perl -lane 'print $F[0]') -gt 0 ]]; then
    last_test=$(tail -n1 $csv | cut -d',' -f1)
    for test in $raw_tests; do
        if [[ $test > $raw_tests ]]; then
            tests+=" $test"
        else 
            echo "skipping $test"
        fi

    done
else
    tests=$raw_tests
fi


function recorded_run() {
    local func="$1"
    local testname="$2"
    if [[ $func == compile ]]; then
        echo "======= compiling $TEST ======="
        timeout --signal=INT --preserve-status $timelimit make $TEST
    elif [[ $func == run ]]; then
        exe=$(fd $testname)
        if [[ -z $exe ]]; then
            return 0; # possibly a compile test
        fi
        echo "======= running ./$exe ======="
        timeout --signal=INT --preserve-status $timelimit ./$exe
    else
        echo "ERROR: unknown func=$func"
        exit 1
    fi
    res=$?
    if [[ $res -eq 124 ]]; then
        echo "======= $func timed out ======="
        if [[ $func == compile ]]; then
            compile_timed_out=1
        else
            exec_timed_out=1
        fi
    elif [[ $res -ne 0 ]]; then
        if [[ $func == compile ]]; then
            compile_failed=1
        else
            exec_failed=1
        fi
    fi
    return $res
}

for test in $tests; do
    compile_timed_out=0
    exec_timed_out=0
    compile_failed=0
    exec_failed=0
    echo "======= starting $test ======="
    if [[ $res -eq 130 ]]; then
        echo "Caught SIGINT during compile"
        exit 130
    elif [[ $res -eq 0 ]]; then 
        recorded_run run $test
        res=$?
        if [[ $res -eq 130 ]]; then
            echo "Caught SIGINT during run"
            exit 130
        fi
    fi
    echo "$test, $arch, $compile_timed_out, $compile_failed, $exec_timed_out, $exec_failed" >> $csv
done
