"""Compile the retained density module against an existing GPU Debug build.

Usage: python3 reproduce.py /absolute/path/to/build-Debug [nvfortran]
The build must use nvfortran and include ScaLAPACK and cuSOLVERMp.
The source uses valid host association (Fortran 2008, 16.5.1.4) to pass
its parent's info dummy to a use-associated procedure from purify.
"""
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

build = Path(sys.argv[1]).resolve()
flags_file = (build / 'src/solvers/CMakeFiles/mspoly_solvers.dir/flags.make')
lines = flags_file.read_text().splitlines()
compiler = sys.argv[2] if len(sys.argv) > 2 else next(
    line.split(' with ', 1)[1] for line in lines
    if line.startswith('# compile Fortran with '))
flags = []
for key in ['Fortran_INCLUDES', 'Fortran_FLAGS']:
    flags += shlex.split(next(line.split(' = ', 1)[1] for line in lines
                              if line.startswith(key + ' = ')))
index = flags.index('-module')
flags[index:index + 2] = []
source = Path(__file__).with_name('msp_density_mod.f90').resolve()
with tempfile.TemporaryDirectory(prefix='msp-density-debug-') as dest:
    result = subprocess.run([compiler, *flags, '-module', dest, '-c',
                             str(source), '-o', str(Path(dest) / 'density.o')],
                            cwd=dest, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    print(result.stdout, end='')
    if result.returncode:
        raise SystemExit(result.returncode)
    if 'invalid debug info' in result.stdout:
        print('RESULT: invalid debug information emitted')
        raise SystemExit(1)
print('RESULT: all passed')
