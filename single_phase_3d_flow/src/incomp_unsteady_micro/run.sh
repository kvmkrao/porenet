
#python 3d_openpnm_simple_gen.py 
make clean
rm -rf CMakeFiles/ CMakeCache.txt  cmake_install.cmake
cmake . 
make 
./sinphase  -ksp_type fgmres  -pc_type hypre -ksp_max_it=60000 -ksp_atol=1e-10
#./sinphase  -ksp_type fgmres  -pc_type hypre  -mat_view -ksp_max_it=60000 -ksp_atol=1e-10
#display mass.png 
