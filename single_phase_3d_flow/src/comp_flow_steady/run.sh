./teton  -ksp_type fgmres  -pc_type hypre  -mat_view -ksp_max_it=60000 -ksp_atol=1e-8
./teton  -ksp_type fgmres  -pc_type jacobi -mat_view -ksp_max_it=60000 -ksp_atol=1e-8
#./teton  -ksp_type gmres  -pc_type jacobi  -mat_view -ksp_max_it=20000 -ksp_atol=1e-6
#./teton  -ksp_type cg -pc_type hypre -mat_view 
