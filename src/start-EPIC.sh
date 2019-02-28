#!/usr/bin/env bash


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "Type in name of network interface, followed by [ENTER]:(leave empty for default: en0)"

read neti

if ["$neti" == ""]; then
        echo "No network specified using en0 as default network interface"
        neti="en0";
fi
echo "Using network interface: $neti"

docker run --user="jovyan" --add-host="localhost:$(ifconfig $neti | grep inet | grep -v inet6 | awk '{print $2}')" -it --rm -p 8888:8888 -v "$DIR:/home/jovyan/work/input" baderlab/bio-epic start-notebook.sh --NotebookApp.iopub_data_rate_limit='100000000' --NotebookApp.token=''
