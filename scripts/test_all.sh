#!/bin/bash

if ! which timeout 2>&1>/dev/null; then
    echo "ERROR: missing required command timeout (brew install coreutils)"
    exit 1
fi

if ! which fd 2>&1>/dev/null; then
    echo "ERROR: missing required command fd (brew install fd)"
    exit 1
fi

filter=""
dryrun=""
timelimit="5m"
only_compile=""
abort_for_compilation_failure=""

if [[ $(uname) == Darwin ]]; then
    arch=mac
elif [[ $(uname) == Linux ]]; then
    arch=linux
else
    arch=unknown
fi

csv=$arch.csv

function show_help() {
    echo "  -h              this message"
    echo "  -d              dry run"
    echo "  -f FILTER       filter test names"
    echo "  -t LIMIT        limit compilation and runtime [$timelimit]"
    echo "  -o FILE         output file [$csv]"
    echo "  -c              only compile"
    echo "  -C              only compile, stop after the first compilation failure"
}

while getopts "h?df:t:o:ocC" opt; do
  case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    d)
        dryrun=1
        ;;
    f)  
        filter=$OPTARG
        ;;
    t)
        timelimit=$OPTARG
        ;;
    o)
        csv=$OPTARG
        ;;
    c)
        only_compile=1
        ;;
    C)
        only_compile=1
        abort_for_compilation_failure=1
        ;;
  esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift



function list_make_targets() {
    make -qp \
        | awk -F':' '/^[a-zA-Z0-9][^$#\/\t=]*:([^=]|$)/ {split($1,A,/ /);for(i in A)print A[i]}' \
        | sort -u
}

raw_tests=$(list_make_targets | grep "$filter")

# drop tests which we've logged already to the csv file (allowing resumption)
touch $csv
if [[ $(wc -l "$csv" | perl -lane 'print $F[0]') -gt 0 ]]; then
    for test in $raw_tests; do
        if [[ ! $(grep $test $csv) ]]; then
            # echo "adding $test"
            tests+=" $test"
        else 
            echo "skipping $test (logged already in $csv)"
        fi
    done
else
    tests=$raw_tests
fi


function recorded_run() {
    local func="$1"
    local testname="$2"
    if [[ $func == compile ]]; then
        echo "======= compiling $testname ======="
        timeout $timelimit make $testname
    elif [[ $func == run ]]; then
        exe=$(fd $testname)
        if [[ -z $exe ]]; then
            return 0; # possibly a compile test
        fi
        echo "======= running ./$exe ======="
        timeout $timelimit ./$exe
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

    [[ -n $dryrun ]] && continue # if dryrun, skip to next iteraton

    recorded_run compile $test

    [[ -n $abort_for_compilation_failure && $compile_failed == 1 ]] && exit 1

    if [[ -n $only_compile ]]; then
        echo "$test, $arch, $compile_timed_out, $compile_failed, , " >> $csv
    else
        recorded_run run $test
        echo "$test, $arch, $compile_timed_out, $compile_failed, $exec_timed_out, $exec_failed" >> $csv
    fi
done
