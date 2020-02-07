#!/bin/bash
#SBATCH --job-name=daskJob_maxm_n8_p8_day
#SBATCH --partition=high_mem
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --mem=MaxMemPerNode
#SBATCH --qos=normal+
#SBATCH --output=z_slurm-%x-%j-%u.out
#SBATCH --error=z_slurm-%x-%j-%u.out
#SBATCH --wait-all-nodes=1

current_time=$(date "+%Y.%m.%d-%H.%M.%S")
which python
echo "Current Time : $current_time"
master=$(echo $SLURMD_NODENAME)
echo "master: $master"
master_ip=$(getent hosts $master  | awk '{ print $1 }')
echo "master_ip: $master_ip"
echo "delete old schedulers"
ps aux | grep -ie dask-scheduler | awk '{print $2}' | xargs kill -9
sleep 3
echo "start new scheduler"
dask-scheduler &
sleep 30

#srun python modis.py $master_ip 

nodes=$(echo $SLURM_NODELIST | cut -d'[' -f 2 | cut -d']' -f 1)
echo "nodes: $nodes"
echo $(echo $nodes| cut -d'-' -f 2)

if [[ $nodes == *","* ]]; then
  for element in ${nodes//,/ } ; do
    echo "element:$element"
    if [[ $element == *"-"* ]]; then
      start_node=$(echo $element| cut -d'-' -f 1)
      end_node=$(echo $element| cut -d'-' -f 2)
      echo "start_node:$start_node, end_node:$end_node"
      for sub_element in $(seq -f "%03g" $start_node $end_node) ; do
        if [ "cnode$sub_element" != "$master" ]; then
          echo 'branch 1'
          worker_node=$(echo 'cnode'$sub_element)
          echo $worker_node
          echo 'cnode'$sub_element
          sleep 5
          script="/umbc/xfs1/cybertrn/common/Softwares/anaconda3/bin/dask-worker $master_ip:8786"
          echo "worker script: ssh ${worker_node} ${script}"
          ssh ${worker_node} /umbc/xfs1/jianwu/common/MODIS_Aggregation/jianwu-code/github-dask/MODIS_Aggregation/benchmarking/dask-environment-setup/start_worker.sh ${master_ip} &
          sleep 5
          echo "exit worker node: $worker_node"
        fi
      done
    else
      if [ "cnode$element" != "$master" ]; then
        echo 'branch 2'
        worker_node=$(echo 'cnode'$element)
        echo $worker_node
        echo 'cnode'$element
        sleep 5
        script="/umbc/xfs1/cybertrn/common/Softwares/anaconda3/bin/dask-worker $master_ip:8786"
        echo "worker script: ssh ${worker_node} ${script}"
        ssh ${worker_node} /umbc/xfs1/jianwu/common/MODIS_Aggregation/jianwu-code/github-dask/MODIS_Aggregation/benchmarking/dask-environment-setup/start_worker.sh ${master_ip} &
        sleep 5
        echo "exit worker node: $worker_node"
      fi
    fi
  done
else
  start_node=$(echo $nodes| cut -d'-' -f 1)
  end_node=$(echo $nodes| cut -d'-' -f 2)
  echo "start_node:$start_node, end_node:$end_node"
  for sub_element in $(seq -f "%03g" $start_node $end_node) ; do
    if [ "cnode$sub_element" != "$master" ]; then
          echo 'branch 3'
          worker_node=$(echo 'cnode'$sub_element)
          echo $worker_node
          echo "ssh into worker node: $worker_node"
          sleep 5
          script="/umbc/xfs1/cybertrn/common/Softwares/anaconda3/bin/dask-worker $master_ip:8786"
          echo "worker script: ssh ${worker_node} ${script}"
          ssh ${worker_node} /umbc/xfs1/jianwu/common/MODIS_Aggregation/jianwu-code/github-dask/MODIS_Aggregation/benchmarking/dask-environment-setup/start_worker.sh ${master_ip} &
          sleep 5
          echo "exit worker node: $worker_node"
	  #dask-worker $worker_node:8786
    fi
  done
fi

which python
sleep 20
echo "Starting running python file"
/umbc/xfs1/cybertrn/common/Softwares/anaconda3/bin/python /home/jianwu/jianwu_common/MODIS_Aggregation/jianwu-code/github-dask/MODIS_Aggregation/benchmarking/xarray-dask/monthly-aggregation-file-level-dask.py $master_ip 
echo "call kill worker script"

if [[ $nodes == *","* ]]; then
  for element in ${nodes//,/ } ; do
    echo "element:$element"
    if [[ $element == *"-"* ]]; then
      start_node=$(echo $element| cut -d'-' -f 1)
      end_node=$(echo $element| cut -d'-' -f 2)
      echo "start_node:$start_node, end_node:$end_node"
      for sub_element in $(seq -f "%03g" $start_node $end_node) ; do
        if [ "cnode$sub_element" != "$master" ]; then
          echo 'branch 1'
          worker_node=$(echo 'cnode'$sub_element)
          echo $worker_node
          echo 'cnode'$sub_element
          sleep 5
          echo "worker script: ssh ${worker_node}"
          ssh ${worker_node} /umbc/xfs1/jianwu/common/MODIS_Aggregation/jianwu-code/github-dask/MODIS_Aggregation/benchmarking/dask-environment-setup/kill_worker.sh ${master_ip} &
          sleep 5
          echo "exit worker node: $worker_node"
        fi
      done
    else
      if [ "cnode$element" != "$master" ]; then
        echo 'branch 2'
        worker_node=$(echo 'cnode'$element)
        echo $worker_node
        echo 'cnode'$element
        sleep 5
        echo "worker script: ssh ${worker_node}"
        ssh ${worker_node} /umbc/xfs1/jianwu/common/MODIS_Aggregation/jianwu-code/github-dask/MODIS_Aggregation/benchmarking/dask-environment-setup/kill_worker.sh ${master_ip} &
        sleep 5
        echo "exit worker node: $worker_node"
      fi
    fi
  done
else
  start_node=$(echo $nodes| cut -d'-' -f 1)
  end_node=$(echo $nodes| cut -d'-' -f 2)
  echo "start_node:$start_node, end_node:$end_node"
  for sub_element in $(seq -f "%03g" $start_node $end_node) ; do
    if [ "cnode$sub_element" != "$master" ]; then
          echo 'branch 3'
          worker_node=$(echo 'cnode'$sub_element)
          echo $worker_node
          echo "ssh into worker node: $worker_node"
          sleep 5
          echo "worker script: ssh ${worker_node}"
          ssh ${worker_node} /umbc/xfs1/jianwu/common/MODIS_Aggregation/jianwu-code/github-dask/MODIS_Aggregation/benchmarking/dask-environment-setup/kill_worker.sh ${master_ip} &
          sleep 5
          echo "exit worker node: $worker_node"
          #dask-worker $worker_node:8786
    fi
  done
fi

echo "call kill_worker finishes"
echo "End running python file"
