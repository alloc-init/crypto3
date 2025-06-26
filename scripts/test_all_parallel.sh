#!/bin/bash

set -euo pipefail

timelimit="5m"
max_jobs=$(nproc)  # number of parallel jobs

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

for cmd in timeout fd flock nproc; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: missing required command $cmd"
        exit 1
    fi
done

list_make_targets() {
    make -qp \
        | awk -F':' '/^[a-zA-Z0-9][^$#\/\t=]*:([^=]|$)/ {split($1,A,/ /);for(i in A)print A[i]}' \
        | sort -u
}

# Get all test targets
mapfile -t raw_tests < <(list_make_targets | grep _test | sort | sed 's/block_aria_test//')

# If csv is empty, create it and write the header
if [[ ! -s "$csv" ]]; then
    echo "test,arch,compile_timed_out,compile_failed,exec_timed_out,exec_failed" > "$csv"
    tests=("${raw_tests[@]}")
else
    last_test=$(tail -n1 "$csv" | cut -d',' -f1)
    tests=()
    for test in "${raw_tests[@]}"; do
        if [[ $test > $last_test ]]; then
            tests+=("$test")
        else
            echo "skipping $test"
        fi
    done
fi

run_test() {
    local test="$1"
    local compile_timed_out=0
    local exec_timed_out=0
    local compile_failed=0
    local exec_failed=0

    echo "======= starting $test ======="

    timeout "$timelimit" make "$test"
    res=$?
    if [[ $res -eq 124 ]]; then
        compile_timed_out=1
        echo "Compile timed out for $test"
    elif [[ $res -ne 0 ]]; then
        compile_failed=1
        echo "Compile failed for $test"
    fi

    if [[ $res -eq 0 ]]; then
        exe=$(fd "$test" . | head -n1)
        if [[ -n $exe ]]; then
            timeout "$timelimit" ./"$exe"
            res=$?
            if [[ $res -eq 124 ]]; then
                exec_timed_out=1
                echo "Run timed out for $test"
            elif [[ $res -ne 0 ]]; then
                exec_failed=1
                echo "Run failed for $test"
            fi
        fi
    fi

    # flock protects csv write from concurrent access
    flock "$csv.lock" -c "echo \"$test,$arch,$compile_timed_out,$compile_failed,$exec_timed_out,$exec_failed\" >> \"$csv\""
}

pids=()
declare -A pid_to_test

for test in "${tests[@]}"; do
    run_test "$test" &
    pid=$!
    pids+=("$pid")
    pid_to_test[$pid]="$test"
    while (( ${#pids[@]} >= max_jobs )); do
        for i in "${!pids[@]}"; do
            if ! kill -0 "${pids[$i]}" 2>/dev/null; then
                wait "${pids[$i]}"
                status=$?
                testname=${pid_to_test[${pids[$i]}]}
                if [[ $status -eq 0 ]]; then
                    echo "[$testname finished]"
                else
                    echo "[$testname failed, exit code $status]"
                fi
                unset 'pids[i]'
                unset 'pid_to_test[${pids[$i]}]'
            fi
        done
        sleep 0.1
    done
done

for pid in "${pids[@]}"; do
    wait "$pid"
    status=$?
    testname=${pid_to_test[$pid]}
    if [[ $status -eq 0 ]]; then
        echo "[$testname finished]"
    else
        echo "[$testname failed, exit code $status]"
    fi
done

echo "All tests completed."
