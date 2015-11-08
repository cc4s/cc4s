# @ job_type=MPICH
# @ output = $(job_name).$(jobid).out
# @ error =  $(job_name).$(jobid).err
# @ class = 72h
# @ restart = no
# @ node = 2
# @ total_tasks = 32
# @ queue
#mpiexec.hydra -bootstrap ll ~/Projects/DistributedCoupledCluster/cc4s/bin/Cc4s -no 8 -nv 32 -nG 1
#mpiexec.hydra -bootstrap ll ~/Projects/DistributedCoupledCluster/ctf/bin/ccsd -no 8 -nv 32 -niter 2
#mpiexec.hydra -bootstrap ll ~/Projects/DistributedCoupledCluster/cc4s/bin/Cc4s -no 50 -nv 700 -nG 1000
mpiexec.hydra -bootstrap ll ~/Projects/DistributedCoupledCluster/cc4s/bin/Cc4s -storeV 0 -niter 100
