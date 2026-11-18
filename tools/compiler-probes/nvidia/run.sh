#!/bin/bash
# Build and run every r*.f90 reproducer in this directory with each given
# compiler, at -O0, -O2 and the compiler's standard-conformance mode.
# Usage: run.sh [compiler...]   (default: nvfortran gfortran)
# Writes one line per reproducer and flag set to results.txt.
cd "$(dirname "$0")"
compilers=${*:-nvfortran gfortran}
out=results.txt; : > $out
log() { echo "$*" | tee -a $out; }
for cc in $compilers; do
  version=$($cc --version 2>&1 | grep -m1 -E 'nvfortran|GNU Fortran')
  if [[ $version == *nvfortran* ]]; then
    flag_sets=("-O0" "-O2" "-O0 -Mstandard")
  else
    flag_sets=("-O0 -std=f2008" "-O2 -std=f2008")
  fi
  for flags in "${flag_sets[@]}"; do
    d=$(mktemp -d)
    log "#### $(basename $cc) [$flags] $version"
    for f in r*.f90; do
      p=${f%.f90}
      if (cd $d && $cc $flags $OLDPWD/$f -o $p > $p.log 2>&1); then
        result=$( (cd $d && timeout 20 ./$p 2>&1); echo "exit=$?")
        log "  $p: $(echo "$result" | grep -E 'RESULT|exit=' | tr '\n' ' ')"
      else
        log "  $p: COMPILE ERROR: $(grep -m1 -E 'S-[0-9]|Error' $d/$p.log | cut -c1-120)"
      fi
    done
    rm -rf $d
  done
done
