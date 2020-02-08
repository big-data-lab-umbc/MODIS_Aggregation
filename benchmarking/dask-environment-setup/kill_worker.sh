echo "start kill_worker"
echo "scheduler node: $1"
echo "Delete old worker"
ps aux | grep -ie dask-worker | awk '{print $2}' | xargs kill -9
hostname
date
echo "finish kill_worker"
