#!/bin/bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 --debug | --release" >&2
    exit 2
fi

case "$1" in
    --debug)
        preset="debug"
        ;;
    --release)
        preset="release"
        ;;
    *)
        echo "Usage: $0 --debug | --release" >&2
        exit 2
        ;;
esac

script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_directory}"

cmake --preset "${preset}"
cmake --build --preset "${preset}" --parallel 2

"${script_directory}/build/${preset}/Deterministic_FVM_Code"
