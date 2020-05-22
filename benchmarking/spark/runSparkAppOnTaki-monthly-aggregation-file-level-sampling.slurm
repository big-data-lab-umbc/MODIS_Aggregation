#!/bin/bash
#SBATCH --partition=high_mem
#SBATCH --exclusive
#SBATCH --qos=normal+
#SBATCH --output=slurm-%x-%j-%u.out
#SBATCH --error=slurm-%x-%j-%u.out

#Demo command to run the same slurm on 2 nodes with sampling rate as 3: node=2 && sampling=3 && sbatch --export=sampling=$sampling,node=$node --job-name=$node-node-$sampling-sampling --nodes=$node runSparkAppOnTaki-monthly-aggregation-file-level-sampling.slurm

export CONDA_PREFIX=/umbc/xfs1/cybertrn/common/Softwares/anaconda3/
export PYSPARK_PYTHON=$CONDA_PREFIX/bin/python
export PYSPARK_DRIVER_PYTHON=$CONDA_PREFIX/bin/python

SPARK=/umbc/xfs1/cybertrn/common/Softwares/spark/spark-2.4.0-bin-hadoop2.7
MY_SPARK=/home/jianwu/jianwu_common/pei_spark/my-spark-jianwu
SPARK_PYTHON_FILE=/home/jianwu/jianwu_common/MODIS_Aggregation/jianwu-code/github/MODIS_Aggregation/benchmarking/spark/monthly-aggregation-file-level-spark-sampling.py $sampling $SLURM_NNODES


current_time=$(date "+%Y.%m.%d-%H.%M.%S")
echo "Current Time : $current_time"
MY_SPARK=$MY_SPARK/spark-slurm-job-$SLURM_JOB_ID-$current_time
mkdir -p $MY_SPARK
EXE_LOG_PATH_MY_SPARK=$MY_SPARK/logs
cp -r $SPARK/conf $MY_SPARK
export SPARK_CONF_DIR=$MY_SPARK/conf
mkdir -p $EXE_LOG_PATH_MY_SPARK


#Step 2: Update slaves file at $SPARK/conf based on the nodes allocated from job scheduler.

echo "node number: $SLURM_NNODES"
cat /dev/null > $MY_SPARK/conf/slaves
master=$(echo $SLURMD_NODENAME)

nodes=$(echo $SLURM_NODELIST | cut -d'[' -f 2 | cut -d']' -f 1)
echo "nodes: $nodes"
echo $(echo $nodes| cut -d'-' -f 2)
if [[ $nodes == *","* ]]; then
  for element in ${nodes//,/ } ; do
    echo "element:$element"
    if [[ $element == *"-"* ]]; then
      start_node=$(echo $element| cut -d'-' -f 1)
      end_node=$(echo $element| cut -d'-' -f 2)
      #echo "start_node:$start_node, end_node:$end_node"
      for sub_element in $(seq -f "%03g" $start_node $end_node) ; do
        # add master node to slaves as well
        #if [ "cnode$sub_element" != "$master" ]; then
        echo 'cnode'$sub_element >> $MY_SPARK/conf/slaves
        #fi
      done
    else
      # add master node to slaves as well
      #if [ "cnode$element" != "$master" ]; then
      echo 'cnode'$element >> $MY_SPARK/conf/slaves
      #fi
    fi
  done
else
  start_node=$(echo $nodes| cut -d'-' -f 1)
  end_node=$(echo $nodes| cut -d'-' -f 2)
  for sub_element in $(seq -f "%03g" $start_node $end_node) ; do
    # add master node to slaves as well
    #if [ "cnode$sub_element" != "$master" ]; then
    echo 'cnode'$sub_element >> $MY_SPARK/conf/slaves
    #fi
  done
fi  

echo "slaves: $(cat $MY_SPARK/conf/slaves)"

echo $(egrep --color 'Mem|Cache|Swap' /proc/meminfo)
echo $(ulimit -a)



#Step 3: Start/deploy Spark on all nodes allocated
$SPARK/sbin/stop-all.sh
sleep 5

ulimit -c unlimited
#$SPARK/sbin/start-master.sh
$SPARK/sbin/start-all.sh
sleep 5

host=$(hostname)
export _JAVA_OPTIONS="-Xmx20g"

#Step 4: Submit your Spark job and wait for its finish
echo "spark submit command: $SPARK/bin/spark-submit --total-executor-cores $((32*$SLURM_NNODES)) --executor-cores=2 --master spark://$master:7077 --driver-memory 20g --executor-memory 20g --conf spark.driver.maxResultSize=0 $SPARK_PYTHON_FILE"

echo "log files of this job can be found at $MY_SPARK"

echo "time $SPARK/bin/spark-submit --total-executor-cores $((32*$SLURM_NNODES)) --executor-cores=2 --master spark://$master:7077 --driver-memory 20g --executor-memory 20g --conf spark.driver.maxResultSize=0 $SPARK_PYTHON_FILE" > $EXE_LOG_PATH_MY_SPARK/$current_time-script-log.txt

(time $SPARK/bin/spark-submit --total-executor-cores $((32*$SLURM_NNODES)) --executor-cores=2 --master spark://$master:7077 --driver-memory 20g --executor-memory 20g --conf spark.driver.maxResultSize=0 $SPARK_PYTHON_FILE > $EXE_LOG_PATH_MY_SPARK/$current_time-log.txt) 2> $EXE_LOG_PATH_MY_SPARK/$current_time-time.txt

#Step 5: Stop Spark at the end
sleep 5
$SPARK/sbin/stop-all.sh
