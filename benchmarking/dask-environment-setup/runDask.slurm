#!/bin/bash
#SBATCH --job-name=daskJob_maxm_n8_p16_day
#SBATCH --partition=high_mem
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --mem=MaxMemPerNode
#SBATCH --qos=normal+
#SBATCH --output=dask_slurm-%x-%j-%u.out
#SBATCH --error=dask_slurm-%x-%j-%u.out
#SBATCH --wait-all-nodes=1


export CONDA_PREFIX=/umbc/xfs1/cybertrn/common/Softwares/anaconda3/

DASK_SLURM_PATH=/umbc/xfs1/jianwu/common/MODIS_Aggregation/jianwu-code/github-dask/MODIS_Aggregation/benchmarking/dask-environment-setup/
DASK_PYTHON_FILE=/umbc/xfs1/jianwu/common/MODIS_Aggregation/jianwu-code/github/MODIS_Aggregation/benchmarking/dask-environment-setup/modis-file-level.py

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

echo "delete old worker on master before run py file"
ps aux | grep -ie dask-worker | awk '{print $2}' | xargs kill -9
echo "start new worker on master"
ssh ${master} $DASK_SLURM_PATH/start_worker.sh ${master_ip} &
sleep 5

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
          script="$CONDA_PREFIX/bin/dask-worker $master_ip:8786"
          echo "worker script: ssh ${worker_node} ${script}"
          ssh ${worker_node} $DASK_SLURM_PATH/start_worker.sh ${master_ip} &
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
        script="$CONDA_PREFIX/bin/dask-worker $master_ip:8786"
        echo "worker script: ssh ${worker_node} ${script}"
        ssh ${worker_node} $DASK_SLURM_PATH/start_worker.sh ${master_ip} &
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
          script="$CONDA_PREFIX/bin/dask-worker $master_ip:8786"
          echo "worker script: ssh ${worker_node} ${script}"
          ssh ${worker_node} $DASK_SLURM_PATH/start_worker.sh ${master_ip} &
          sleep 5
          echo "exit worker node: $worker_node"
	  #dask-worker $worker_node:8786
    fi
  done
fi

which python
sleep 20
echo "Starting running python file"
$CONDA_PREFIX/bin/python $DASK_PYTHON_FILE/modis-file-level.py $master_ip 
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
          ssh ${worker_node} $DASK_SLURM_PATH/kill_worker.sh ${master_ip} &
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
        ssh ${worker_node} $DASK_SLURM_PATH/kill_worker.sh ${master_ip} &
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
          ssh ${worker_node} $DASK_SLURM_PATH/kill_worker.sh ${master_ip} &
          sleep 5
          echo "exit worker node: $worker_node"
          #dask-worker $worker_node:8786
    fi
  done
fi

echo "delete old schedulers after py run"
ps aux | grep -ie dask-scheduler | awk '{print $2}' | xargs kill -9

echo "delete old worker on master after run py file"
ps aux | grep -ie dask-worker | awk '{print $2}' | xargs kill -9

echo "call kill_worker finishes"
echo "End running python file"
