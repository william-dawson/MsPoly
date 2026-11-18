#!/bin/bash
# Line coverage of the library sources by the test suite (gfortran + gcov).
# Usage: tools/coverage.sh [build-dir [cmake-args...]]   (default: build-coverage)
set -e
root=$(cd "$(dirname "$0")/.." && pwd)
dir=${1:-$root/build-coverage}
shift || true
gcov=${GCOV:-$(command -v gcov-16 || command -v gcov)}
cmake -S "$root" -B "$dir" -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_Fortran_FLAGS="--coverage -O0" -DCMAKE_EXE_LINKER_FLAGS="--coverage" "$@" > /dev/null
cmake --build "$dir" > /dev/null
find "$dir" \( -name '*.gcda' -o -name '*.gcov' \) -delete
ctest --test-dir "$dir" > /dev/null || { echo "tests failed"; exit 1; }
cd "$dir"
for obj in $(find src -name '*.gcda'); do
  (cd "$(dirname "$obj")" && "$gcov" "$(basename "$obj")" > /dev/null)
done
echo "source                          lines   covered"
for f in $(find src -name '*.gcov' | grep -E 'msp_.*\.(f90|fypp)\.gcov$' | sort -u); do
  total=$(grep -c -E '^ +([0-9]+\*?|#####):' "$f" || true)
  missed=$(grep -c -E '^ +#####:' "$f" || true)
  printf "%-30s %6d %8.1f%%\n" $(basename "$f" .gcov) $total \
    $(echo "scale=1; 100*($total-$missed)/$total" | bc)
done
echo
echo "Unexecuted statements (file:line), ignoring lines of continued statements,"
echo "where gcov attributes counts to the first line only:"
for f in $(find src -name '*.gcov' | grep -E 'msp_.*\.(f90|fypp)\.gcov$' | sort -u); do
  name=$(basename "$f" .gcov)
  source=$(find "$root/src" -name "$name")
  awk -F: '/^ +#####:/ {gsub(/ /, "", $2); print $2}' "$f" | sort -n -u | while read n; do
    line=$(sed -n "${n}p" "$source")
    previous=$(sed -n "$((n - 1))p" "$source")
    case "$line$previous" in
      *"&"*|*"end module"*) ;;
      *) echo "  $name:$n: $(echo $line | cut -c1-80)" ;;
    esac
  done
done
