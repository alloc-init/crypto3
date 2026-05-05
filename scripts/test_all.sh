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
parallel_jobs="${TEST_ALL_JOBS:-1}"

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
    echo "  -j JOBS         tests to run in parallel [$parallel_jobs]"
}

while getopts "h?df:t:o:ocsj:" opt; do
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
    j)
        parallel_jobs=$OPTARG
        ;;
  esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

if ! [[ $parallel_jobs =~ ^[0-9]+$ ]] || [[ $parallel_jobs -lt 1 ]]; then
    echo "ERROR: -j requires a positive integer"
    exit 1
fi

if [[ -n $abort_on_failure && $parallel_jobs -gt 1 ]]; then
    echo "ERROR: -s is only supported with -j 1"
    exit 1
fi



function list_make_targets() {
    make -qp \
        | awk -F':' '/^[a-zA-Z0-9][^$#\/\t=]*:([^=]|$)/ {split($1,A,/ /);for(i in A)print A[i]}' \
        | grep -e "_test$" \
        | sort -u
}

raw_tests=$(list_make_targets | grep "$filter")

function csv_has_test() {
    local testname="$1"
    awk -F', *' -v testname="$testname" '
        $1 == testname { found = 1; exit }
        END { exit(found ? 0 : 1) }
    ' "$csv"
}

function csv_has_runtime_result() {
    local testname="$1"
    awk -F', *' -v testname="$testname" '
        $1 == testname { found = 1; runtime_result = $6; exit }
        END { exit(found && runtime_result != "" ? 0 : 1) }
    ' "$csv"
}

# drop tests which we've logged already to the csv file (allowing resumption)
touch "$csv"
if [[ $(wc -l "$csv" | perl -lane 'print $F[0]') -gt 0 ]]; then
    for test in $raw_tests; do
        if ! csv_has_test "$test"; then
            # echo "adding $test"
            tests+=" $test"
        elif [[ -z $only_compile ]] && ! csv_has_runtime_result "$test"; then
            # if we are running tests and this test has been compiled but not tested, include it in the list
            tests+=" $test"
        else 
            echo "skipping $test (logged already in $csv)"
        fi
    done
else
    tests=$raw_tests
fi

if [[ -z $tests ]]; then
    echo "no tests selected"
    exit 0
fi

if [[ -n $dryrun ]]; then
    for test in $tests; do
        echo "======= starting $test ======="
    done
    exit 0
fi

results_dir=$(mktemp -d "${TMPDIR:-/tmp}/test_all.XXXXXX")
rows_dir="$results_dir/rows"
logs_dir="$csv.logs"
mkdir -p "$rows_dir" "$logs_dir"

function stop_running_jobs() {
    local pid
    for pid in $(jobs -rp); do
        kill "$pid" 2>/dev/null
    done
}

function cleanup() {
    if [[ -n $results_dir && -d $results_dir ]]; then
        rm -rf "$results_dir"
    fi
}
trap cleanup EXIT
trap 'stop_running_jobs; cleanup; exit 130' INT
trap 'stop_running_jobs; cleanup; exit 143' TERM

function compile() {
    local testname="$1"
    compile_failed=0
    compile_timed_out=0
    echo "======= compiling $testname ======="
    timeout -f "$timelimit" make "$testname"
    res=$?
    if [[ $res -eq 124 ]]; then
        echo "======= compiling $testname timed out ======="
        compile_timed_out=1
    elif [[ $res -ne 0 ]]; then
        compile_failed=1
    fi
}

function run_test() {
    local testname="$1"
    test_failed=0
    test_timed_out=0
    exe=$(fd "$testname" | head -n 1)
    if [[ -z $exe ]]; then
        # possibly a compile-only test
        return
    fi
    echo "======= running ./$exe ======="
    timeout -f "$timelimit" "./$exe"
    res=$?
    if [[ $res -eq 124 ]]; then
        echo "======= running $testname timed out ======="
        test_timed_out=1
    elif [[ $res -ne 0 ]]; then
        test_failed=1
    fi
}

function record_result_row() {
    local testname="$1"
    printf '%s, %s, %s, %s, %s, %s\n' \
        "$testname" "$arch" "$compile_timed_out" "$compile_failed" "$test_timed_out" "$test_failed" \
        > "$rows_dir/$testname.csv"
}

function run_one() {
    local testname="$1"
    compile_timed_out=""
    compile_failed=""
    test_timed_out=""
    test_failed=""

    echo "======= starting $testname ======="

    compile "$testname"

    if [[ $compile_failed != 1 && $compile_timed_out != 1 && -z $only_compile ]]; then
        run_test "$testname"
    fi

    record_result_row "$testname"

    if [[ $compile_failed == 1 || $compile_timed_out == 1 || $test_failed == 1 || $test_timed_out == 1 ]]; then
        return 1
    fi

    return 0
}

function merge_results() {
    local new_rows="$results_dir/new.csv"
    local tmp_csv

    : > "$new_rows"
    for test in $tests; do
        if [[ -f "$rows_dir/$test.csv" ]]; then
            cat "$rows_dir/$test.csv" >> "$new_rows"
        fi
    done

    if [[ ! -s $new_rows ]]; then
        return
    fi

    tmp_csv=$(mktemp "${csv}.tmp.XXXXXX") || exit 1
    awk -F', *' '
        NR == FNR {
            if (!($1 in row)) {
                order[++count] = $1
            }
            row[$1] = $0
            next
        }

        $1 in row {
            if (!($1 in emitted)) {
                print row[$1]
                emitted[$1] = 1
            }
            next
        }

        { print }

        END {
            for (i = 1; i <= count; ++i) {
                testname = order[i]
                if (!(testname in emitted)) {
                    print row[testname]
                }
            }
        }
    ' "$new_rows" "$csv" > "$tmp_csv" && mv "$tmp_csv" "$csv"
}

function run_serial() {
    local failures=0
    local status
    local log_file

    for test in $tests; do
        log_file="$logs_dir/$test.log"
        run_one "$test" > "$log_file" 2>&1
        status=$?

        if [[ $status -ne 0 ]]; then
            failures=$((failures + 1))
            echo "======= $test failed (log: $log_file) ======="
            if [[ -n $abort_on_failure ]]; then
                merge_results
                exit 1
            fi
        else
            echo "======= $test finished (log: $log_file) ======="
        fi
    done

    if [[ $failures -gt 0 ]]; then
        echo "======= $failures test(s) failed or timed out ======="
    fi
}

function running_jobs_count() {
    jobs -rp | wc -l | tr -d ' '
}

function run_parallel() {
    local pids=()
    local pid_tests=()
    local pid
    local status
    local failures=0
    local log_file
    local i

    for test in $tests; do
        while [[ $(running_jobs_count) -ge $parallel_jobs ]]; do
            sleep 0.2
        done

        log_file="$logs_dir/$test.log"
        echo "======= starting $test (log: $log_file) ======="
        run_one "$test" > "$log_file" 2>&1 &
        pid=$!
        pids+=("$pid")
        pid_tests+=("$test")
    done

    for i in "${!pids[@]}"; do
        status=0
        wait "${pids[$i]}" || status=$?
        if [[ $status -ne 0 ]]; then
            failures=$((failures + 1))
            echo "======= ${pid_tests[$i]} failed (log: $logs_dir/${pid_tests[$i]}.log) ======="
        else
            echo "======= ${pid_tests[$i]} finished (log: $logs_dir/${pid_tests[$i]}.log) ======="
        fi
    done

    if [[ $failures -gt 0 ]]; then
        echo "======= $failures test(s) failed or timed out ======="
    fi
}

if [[ $parallel_jobs -eq 1 ]]; then
    run_serial
else
    run_parallel
fi

merge_results
echo "======= results written to $csv ======="
echo "======= logs written to $logs_dir ======="
