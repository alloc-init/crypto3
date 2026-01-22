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
abort_on_failure=""

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
    echo "  -s              stop after the first failure"
}

while getopts "h?df:t:o:ocs" opt; do
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
    s)
        abort_on_failure=1
        ;;
  esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift



function list_make_targets() {
    make -qp \
        | awk -F':' '/^[a-zA-Z0-9][^$#\/\t=]*:([^=]|$)/ {split($1,A,/ /);for(i in A)print A[i]}' \
        | grep -e "_test$" \
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
        elif [[ -z $only_compile && $(grep $test $csv | perl -F',\s*' -lane 'print $F[5]') == "" ]]; then
            # if we are running tests and this test has been compiled but not tested, include it in the list
            tests+=" $test"
        else 
            echo "skipping $test (logged already in $csv)"
        fi
    done
else
    tests=$raw_tests
fi

function compile() {
    local testname="$1"
    compile_failed=0
    compile_timed_out=0
    echo "======= compiling $testname ======="
    timeout -f $timelimit make $testname
    res=$?
    if [[ $res -eq 124 ]]; then
        echo "======= $func timed out ======="
        compile_timed_out=1
    elif [[ $res -ne 0 ]]; then
        compile_failed=1
    fi
}

function test() {
    local testname="$1"
    test_failed=0
    test_timed_out=0
    exe=$(fd $testname)
    if [[ -z $exe ]]; then
        # possibly a compile-only test
        return
    fi
    echo "======= running ./$exe ======="
    timeout -f $timelimit ./$exe
    res=$?
    if [[ $res -eq 124 ]]; then
        echo "======= $func timed out ======="
        test_timed_out=1
    elif [[ $res -ne 0 ]]; then
        test_failed=1
    fi
}

function record_results() {
    new_entry="$test, $arch, $compile_timed_out, $compile_failed, $test_timed_out, $test_failed"
    old_entry_lineno=$(grep -n $test $csv | cut -d':' -f1)
    if [[ -z $old_entry_lineno ]]; then 
        # no entry yet, add a new one
        echo "$new_entry" >> $csv
    else
        # replace existing line with new line
        sed -i'' -e "${old_entry_lineno}s/.*/$new_entry/" $csv
    fi
}

for test in $tests; do
    compile_timed_out=""
    compile_failed=""
    test_timed_out=""
    test_failed=""

    echo "======= starting $test ======="
    [[ -n $dryrun ]] && continue # if dryrun, skip to next iteraton

    compile $test
    [[ -n $abort_on_failure && $compile_failed == 1 ]] && exit 1

    if [[ ! $only_compile ]]; then
        test $test
        [[ -n $abort_on_failure && $test_failed == 1 ]] && exit 1
    fi

    record_results
done
