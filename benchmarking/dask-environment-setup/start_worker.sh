echo "start start_worker"
echo "scheduler node: $1"
echo "Delete old worker"
ps aux | grep -ie dask-worker | awk '{print $2}' | xargs kill -9
echo "start new worker"
date
hostname
#start 16 parallel processes and 20 GB per process
/umbc/xfs1/cybertrn/common/Softwares/anaconda3/bin/dask-worker $1:8786 --nprocs 16 --nthreads 1 --memory-limit 20.00GB --death-timeout 60
date
hostname
echo "finish start_worker"
