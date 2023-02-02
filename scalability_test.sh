#/bin/bash

meshSize=$1

mpirun --oversubscribe -np 1 DOMAIN_DECOMPOSITION ExportFolder:string="./Run_N1" MeshParameter:double=$meshSize SchurSolverMaxIterations:uint=100 SchurSolverTolerance:double=1e-16 SchurSolverType:bool=1
mpirun --oversubscribe -np 4 DOMAIN_DECOMPOSITION ExportFolder:string="./Run_N4" MeshParameter:double=$meshSize SchurSolverMaxIterations:uint=100 SchurSolverTolerance:double=1e-16 SchurSolverType:bool=1
mpirun --oversubscribe -np 16 DOMAIN_DECOMPOSITION ExportFolder:string="./Run_N16" MeshParameter:double=$meshSize SchurSolverMaxIterations:uint=100 SchurSolverTolerance:double=1e-16 SchurSolverType:bool=1
