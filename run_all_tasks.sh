#!/bin/bash

joinByString() {
  local separator="$1"
  shift
  local first="$1"
  shift
  printf "%s" "$first" "${@/#/$separator}"
}

allTasks=`ls ./tasks/`

# echo $allTasks

taskArray=`joinByString , $allTasks`

sbatch --array=$taskArray ./execute_flow.sh
# sbatch --array=1-1000 ./execute_flow.sh