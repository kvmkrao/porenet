
#python 3d_openpnm_simple_gen.py 
#cmake . 
#make 
./sinphase  -ksp_type fgmres  -pc_type hypre -ksp_max_it=60000 -ksp_atol=1e-20
#./sinphase  -ksp_type fgmres  -pc_type hypre  -mat_view -ksp_max_it=60000 -ksp_atol=1e-10
#display mass.png 
