#!/bin/bash
#
# Setup a batch processing job for telemetered data files from a cabled EK60
# bioacoustic sonar sensor.
#
# C. Wingard 2021-11-10

# Parse the command line inputs, setting the data directories and processing dates
if [ $# -ne 5 ]; then
    echo "$0: required inputs are the site name, the path to the raw and processed data"
    echo "directories, and the starting and ending dates (format is YYYY-MM-DD) of the"
    echo "year to batch process."
    echo ""
    echo "    example: $0 CE04OSPS /home/ooiuser/data/raw/CE04OSPS/PC01B/05-ZPLSCB102/2017 \\ "
    echo "        /home/ooiuser/data/raw/CE04OSPS/PC01B/05-ZPLSCB102/2017/processed \\ "
    echo "        \"2017-01-01\" \"2017-09-16\""
    exit 1
fi
SITE=${1^^}
DATA_DIR=$2
PROC_DIR=$3
START_DATE=`date -u +%Y%m%d -d $4`
END_DATE=`date -u +%Y%m%d -d $5`

# activate the echogram environment and cd to the ooi_zpls_echograms code directory
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate echogram
cd ~/Documents/GitHub/ooi-zpls-echograms/ooi_zpls_echograms

# Set up concurrent parallel processing using 4 cores (equates to 4 weeks)
N=4

# process the data, using 2012-01-01 as the base year for all plots
for d in $(seq $(date -u +%s -d "2012-01-01") +604800 $(date -u +%s -d $END_DATE)); do
    start_date=`date -u +%Y%m%d -d @$d`
    stop_date=`date -u +%Y%m%d -d "$start_date+7days"`
    if [[ $stop_date -gt $START_DATE ]]; then 
        (./zpls_echogram.py -s $SITE -d $DATA_DIR -o $PROC_DIR -dr $start_date $stop_date -zm EK60) &
    fi
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        # there are already $N jobs running in parallel, wait for one to finish
        wait -n
    fi
done
# no more jobs to run, but wait for the last ones to finish
wait
