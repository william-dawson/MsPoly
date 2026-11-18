#!/bin/bash
module load nvhpc-hpcx/26.9 >/dev/null 2>&1
cd $1
echo "== unbound"; mpirun -np 4 --bind-to none ./place 2>&1 | sort
echo "== bound to the GPU's socket"
mpirun -np 4 --bind-to none bash -c 'numactl --cpunodebind=$((OMPI_COMM_WORLD_LOCAL_RANK / 2)) --membind=$((OMPI_COMM_WORLD_LOCAL_RANK / 2)) ./place' 2>&1 | sort
