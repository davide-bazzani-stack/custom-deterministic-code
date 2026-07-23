#!/bin/bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 --debug | --release" >&2
    exit 2
fi

case "$1" in
    --debug)
        preset="test-debug"
        mode="debug"
        ;;
    --release)
        preset="test-release"
        mode="release"
        ;;
    *)
        echo "Usage: $0 --debug | --release" >&2
        exit 2
        ;;
esac

script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repository_root="$(cd -- "${script_directory}/.." && pwd)"
log_directory="${repository_root}/tests/build/${mode}/logs"
timestamp="$(date '+%Y%m%d_%H%M%S')"
stdout_log="${log_directory}/unit_test_${timestamp}.stdout.log"
stderr_log="${log_directory}/unit_test_${timestamp}.stderr.log"

mkdir -p "${log_directory}"

exec > >(tee "${stdout_log}") 2> >(tee "${stderr_log}" >&2)

echo "Standard output log: ${stdout_log}"
echo "Standard error log:  ${stderr_log}"

cd "${repository_root}"

cmake --preset "${preset}"
cmake --build --preset "${preset}" --parallel 2
"${repository_root}/tests/build/${mode}/tests/fvm_unit_tests"
