#!/bin/bash
set -eu
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

check_only="no"
in_place="-i"
output_replacement=""
if [[ "${1:-no}" = "--check" ]]; then
    check_only="yes"
    in_place=""
    output_replacement="--output-replacements-xml"
fi

clang_fmt="/usr/bin/clang-format-13"

files=$(git ls-tree -r --name-only HEAD|grep -E "(\.c$|\.cpp$|\.h$|\.hpp$)"|grep -v "include/vcl")

set +e
# when formatting occures, we see offset="some number"
$clang_fmt $in_place -style=file $output_replacement ${files}|grep "offset"
no_files_formatted=$?
set -e

if [[ no_files_formatted -eq 0 ]]; then
    echo $no_files_formatted
    if [[ "$check_only" = "yes" ]]; then
        echo "Found formatting errors. Run tools/format-repo.sh"

        exit 1
    fi
fi

exit 0
