echo "start start_worker"
echo "scheduler node: $1"
echo "Delete old worker"
ps aux | grep -ie dask-worker | awk '{print $2}' | xargs kill -9
echo "start new worker"
date
hostname
/umbc/xfs1/cybertrn/common/Softwares/anaconda3/bin/dask-worker $1:8786 --nthreads 1 --memory-limit=0 --death-timeout 60
date
hostname
echo "finish start_worker"
