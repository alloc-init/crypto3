#!/bin/bash

if ! which timeout 2>&1>/dev/null; then
    echo "ERROR: missing required command timeout (brew install coreutils)"
    exit 1
fi

timeout_foreground_arg=""
if timeout -f 1 true >/dev/null 2>&1; then
    timeout_foreground_arg="-f"
elif timeout --foreground 1 true >/dev/null 2>&1; then
    timeout_foreground_arg="--foreground"
fi

filter=""
rerun_test=""
dryrun=""
timelimit="10m"
only_compile=""
abort_on_failure=""
print_failure_logs=""
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
    echo "  -r TEST         rerun one exact test and update the csv even if logged"
    echo "  -t LIMIT        limit compilation and runtime [$timelimit]"
    echo "  -o FILE         output file [$csv]"
    echo "  -c              only compile"
    echo "  -s              stop after the first failure"
    echo "  -j JOBS         tests to run in parallel [$parallel_jobs]"
    echo "  -l              print logs for failed tests"
}

while getopts "h?df:r:t:o:ocsj:l" opt; do
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
    r)
        rerun_test=$OPTARG
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
    l)
        print_failure_logs=1
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

if [[ -n $rerun_test && -n $filter ]]; then
    echo "ERROR: -r cannot be combined with -f"
    exit 1
fi

if [[ -z ${NO_COLOR:-} && ( -t 1 || -n ${FORCE_COLOR:-} ) ]]; then
    color_red=$'\033[31m'
    color_green=$'\033[32m'
    color_yellow=$'\033[33m'
    color_blue=$'\033[34m'
    color_bold=$'\033[1m'
    color_reset=$'\033[0m'
else
    color_red=""
    color_green=""
    color_yellow=""
    color_blue=""
    color_bold=""
    color_reset=""
fi



function list_make_targets() {
    make -qp \
        | awk -F':' '/^[a-zA-Z0-9][^$#\/\t=]*:([^=]|$)/ {split($1,A,/ /);for(i in A)print A[i]}' \
        | grep -e "_test$" \
        | sort -u
}

if [[ -n $rerun_test ]]; then
    raw_tests=$(list_make_targets | awk -v testname="$rerun_test" '$0 == testname')
    if [[ -z $raw_tests ]]; then
        echo "ERROR: no test target named $rerun_test"
        exit 1
    fi
else
    raw_tests=$(list_make_targets | grep "$filter")
fi
logs_dir="$csv.logs"

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

function find_test_executable() {
    local testname="$1"
    find . -type f -name "$testname" -print | sed 's#^\./##' | head -n 1
}

function find_test_target_dir() {
    local testname="$1"
    find . -type d -path "*/CMakeFiles/${testname}.dir" -print | sed 's#^\./##' | head -n 1
}

function find_test_build_output() {
    local testname="$1"
    local exe
    local target_dir

    exe=$(find_test_executable "$testname")
    if [[ -n $exe ]]; then
        echo "$exe"
        return
    fi

    target_dir=$(find_test_target_dir "$testname")
    if [[ -n $target_dir ]]; then
        target_dir=${target_dir%/}
        echo "${target_dir%/CMakeFiles/$testname.dir}/$testname"
    fi
}

function test_needs_rebuild() {
    local testname="$1"
    local output
    local target_dir
    local status

    output=$(find_test_build_output "$testname")
    if [[ -z $output ]]; then
        return 1
    fi

    target_dir=$(find_test_target_dir "$testname")
    if [[ -n $target_dir ]]; then
        target_dir=${target_dir%/}
        make -f "$target_dir/build.make" "$target_dir/depend" >/dev/null 2>&1 || return 0
        make -q -f "$target_dir/build.make" "$output" >/dev/null 2>&1
    else
        make -q "$output" >/dev/null 2>&1
    fi

    status=$?

    [[ $status -ne 0 ]]
}

# drop tests which we've logged already to the csv file (allowing resumption)
touch "$csv"
if [[ -n $rerun_test ]]; then
    tests=$raw_tests
elif [[ $(wc -l "$csv" | perl -lane 'print $F[0]') -gt 0 ]]; then
    for test in $raw_tests; do
        if ! csv_has_test "$test"; then
            # echo "adding $test"
            tests+=" $test"
        elif test_needs_rebuild "$test"; then
            echo "rerunning $test (needs rebuild)"
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
    no_tests_selected=1
fi

if [[ -n $dryrun ]]; then
    if [[ -n $no_tests_selected ]]; then
        echo "no tests selected"
        exit 0
    fi

    for test in $tests; do
        echo "======= starting $test ======="
    done
    exit 0
fi

results_dir=$(mktemp -d "${TMPDIR:-/tmp}/test_all.XXXXXX")
rows_dir="$results_dir/rows"
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

function run_with_timeout() {
    local limit="$1"
    shift

    if [[ -n $timeout_foreground_arg ]]; then
        timeout "$timeout_foreground_arg" "$limit" "$@"
    else
        timeout "$limit" "$@"
    fi
}

function compile() {
    local testname="$1"
    compile_failed=0
    compile_timed_out=0
    echo "======= compiling $testname ======="
    run_with_timeout "$timelimit" make "$testname"
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
    exe=$(find_test_executable "$testname")
    if [[ -z $exe ]]; then
        # possibly a compile-only test
        return
    fi
    echo "======= running ./$exe ======="
    run_with_timeout "$timelimit" "./$exe"
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

function join_reasons() {
    local joined=""
    local reason

    for reason in "$@"; do
        if [[ -z $joined ]]; then
            joined="$reason"
        else
            joined="$joined, $reason"
        fi
    done

    echo "$joined"
}

function report_tests_list() {
    if [[ -n $rerun_test ]]; then
        awk -F', *' '$1 != "" && !seen[$1]++ { print $1 }' "$csv"
    else
        printf '%s\n' $raw_tests
    fi
}

function count_report_tests() {
    local count=0
    local test

    for test in $(report_tests_list); do
        count=$((count + 1))
    done

    echo "$count"
}

function print_test_report() {
    local report_count
    local report_label="selected"
    local completed_count=0
    local failed_count=0
    local passed_count=0
    local not_run_count
    local compile_timeout_count=0
    local compile_error_count=0
    local test_timeout_count=0
    local test_error_count=0
    local test
    local row
    local row_testname
    local row_arch
    local row_compile_timed_out
    local row_compile_failed
    local row_test_timed_out
    local row_test_failed
    local reasons
    local reason_color
    local failed_lines=()
    local failed_tests=()

    if [[ -n $rerun_test ]]; then
        report_label="csv rows"
    fi
    report_count=$(count_report_tests)

    for test in $(report_tests_list); do
        if ! csv_has_test "$test"; then
            continue
        fi

        completed_count=$((completed_count + 1))
        row=$(awk -F', *' -v testname="$test" '$1 == testname { print; exit }' "$csv")
        IFS=',' read -r row_testname row_arch row_compile_timed_out row_compile_failed row_test_timed_out row_test_failed <<< "$row"
        row_compile_timed_out=${row_compile_timed_out//[[:space:]]/}
        row_compile_failed=${row_compile_failed//[[:space:]]/}
        row_test_timed_out=${row_test_timed_out//[[:space:]]/}
        row_test_failed=${row_test_failed//[[:space:]]/}

        if [[ $row_compile_timed_out == 1 || $row_compile_failed == 1 || $row_test_timed_out == 1 || $row_test_failed == 1 ]]; then
            local failure_reasons=()
            failed_count=$((failed_count + 1))

            if [[ $row_compile_timed_out == 1 ]]; then
                failure_reasons+=("compilation timeout")
                compile_timeout_count=$((compile_timeout_count + 1))
            fi
            if [[ $row_compile_failed == 1 ]]; then
                failure_reasons+=("compilation error")
                compile_error_count=$((compile_error_count + 1))
            fi
            if [[ $row_test_timed_out == 1 ]]; then
                failure_reasons+=("test timeout")
                test_timeout_count=$((test_timeout_count + 1))
            fi
            if [[ $row_test_failed == 1 ]]; then
                failure_reasons+=("test error")
                test_error_count=$((test_error_count + 1))
            fi

            reasons=$(join_reasons "${failure_reasons[@]}")
            if [[ $reasons == *timeout* ]]; then
                reason_color=$color_yellow
            else
                reason_color=$color_red
            fi
            failed_lines+=("  ${color_red}FAIL${color_reset} $test ${reason_color}[$reasons]${color_reset} log: $logs_dir/$test.log")
            failed_tests+=("$test")
        else
            passed_count=$((passed_count + 1))
        fi
    done

    not_run_count=$((report_count - completed_count))

    echo
    echo "${color_bold}${color_blue}======= test report =======${color_reset}"
    echo "$report_label: $report_count, completed: $completed_count, ${color_green}passed: $passed_count${color_reset}, ${color_red}failed: $failed_count${color_reset}"

    if [[ $not_run_count -gt 0 ]]; then
        echo "${color_yellow}not run: $not_run_count${color_reset}"
    fi

    if [[ $failed_count -eq 0 ]]; then
        echo "${color_green}all completed tests passed${color_reset}"
        return 0
    fi

    echo "${color_bold}failure types:${color_reset} compilation timeouts: $compile_timeout_count, compilation errors: $compile_error_count, test timeouts: $test_timeout_count, test errors: $test_error_count"
    echo "${color_bold}failed tests:${color_reset}"
    printf '%s\n' "${failed_lines[@]}"

    if [[ -n $print_failure_logs ]]; then
        print_failed_test_logs "${failed_tests[@]}"
    fi

    return 1
}

function print_failed_test_logs() {
    local test
    local log_file

    echo
    echo "${color_bold}${color_blue}======= failed test logs =======${color_reset}"

    for test in "$@"; do
        log_file="$logs_dir/$test.log"
        echo
        echo "${color_bold}======= $test ($log_file) =======${color_reset}"

        if [[ -f $log_file ]]; then
            cat "$log_file"
        else
            echo "${color_yellow}missing log: $log_file${color_reset}"
        fi
    done
}

function run_one() {
    local testname="$1"
    compile_timed_out=0
    compile_failed=0
    test_timed_out=0
    test_failed=0

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

if [[ -n $no_tests_selected ]]; then
    echo "no tests selected"
    echo "======= results written to $csv ======="
    echo "======= logs written to $logs_dir ======="
    print_test_report
    exit $?
fi

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
                return 1
            fi
        else
            echo "======= $test finished (log: $log_file) ======="
        fi
    done

    if [[ $failures -gt 0 ]]; then
        echo "======= $failures test(s) failed or timed out ======="
        return 1
    fi

    return 0
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
        return 1
    fi

    return 0
}

if [[ $parallel_jobs -eq 1 ]]; then
    run_serial
    run_status=$?
else
    run_parallel
    run_status=$?
fi

merge_results
echo "======= results written to $csv ======="
echo "======= logs written to $logs_dir ======="
print_test_report
report_status=$?

if [[ $run_status -ne 0 || $report_status -ne 0 ]]; then
    exit 1
fi

exit 0
