#!/usr/bin/env bash
set -euo pipefail

repo_root="${1:-.}"
cd "$repo_root"

json='{"include":['
separator=''

while IFS= read -r test_dir; do
    project_dir="${test_dir%/test}"
    cmake_file="$project_dir/CMakeLists.txt"
    group="${test_dir#libs/}"
    group="${group%/test}"
    group="${group//\//-}"
    project="$(
        sed -nE 's/^[[:space:]]*cm_project\(([A-Za-z0-9_+-]+).*/\1/p' "$cmake_file" \
            | head -n 1
    )"

    if [[ -z $project ]]; then
        echo "ERROR: cannot find cm_project(...) in $cmake_file" >&2
        exit 1
    fi

    json="${json}${separator}{\"group\":\"${group}\",\"project\":\"${project}\",\"build_target\":\"tests-crypto3-${project}\",\"ctest_label\":\"crypto3-${project}\"}"
    separator=','
done < <(find libs -path '*/test/CMakeLists.txt' -print | sed 's#/CMakeLists.txt$##' | sort)

json="${json}]}"

if [[ $json == '{"include":[]}' ]]; then
    echo "ERROR: no test groups found under libs/**/test/CMakeLists.txt" >&2
    exit 1
fi

echo "$json"
