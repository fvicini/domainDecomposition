#/bin/bash

mpirun --oversubscribe -np 1 DOMAIN_DECOMPOSITION ExportFolder:string="./Run_N1" MeshParameter:double=0.0625 SchurSolverMaxIterations:uint=100 SchurSolverTolerance:double=1e-16 SchurSolverType:bool=1
mpirun --oversubscribe -np 4 DOMAIN_DECOMPOSITION ExportFolder:string="./Run_N4" MeshParameter:double=0.0625 SchurSolverMaxIterations:uint=100 SchurSolverTolerance:double=1e-16 SchurSolverType:bool=1
mpirun --oversubscribe -np 16 DOMAIN_DECOMPOSITION ExportFolder:string="./Run_N16" MeshParameter:double=0.0625 SchurSolverMaxIterations:uint=100 SchurSolverTolerance:double=1e-16 SchurSolverType:bool=1
