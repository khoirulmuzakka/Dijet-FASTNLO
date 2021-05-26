Sherpa -f Run_step1_proc.dat >&! Run_step1_proc.log&
rm -rf Results.db Results.db.bak Sherpa_References.tex Analysis.yoda mcgrid
mpiexec -np 4 Sherpa -f Run_step2_inte.dat >&! Run_step2_inte.log&
rm -rf Sherpa_References.tex Analysis*.yoda mcgrid
mpiexec -np 4 Sherpa -f Run_step3_warm.dat > & ! Run_step3_warm.log &
rm -rf Sherpa_References.tex Analysis*.yoda
mpiexec -np 4 Sherpa -f Run_step4_prod.dat -e 16M -R '26763 3746923 10826382 7552378' >&! Run_step4_prod_16M.log
