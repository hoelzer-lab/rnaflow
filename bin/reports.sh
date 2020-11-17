#!/usr/bin/bash

# copy final reports to output dir

DATE=$(date "+%Y-%m-%dT%H%M")

OUTPUT_DIR=$1
REPORTS_DIR=$2

mkdir -p $OUTPUT_DIR/Logs

cp $REPORTS_DIR/execution_report.html $OUTPUT_DIR/Logs/$DATE"_execution_report.html"
cp $REPORTS_DIR/execution_timeline.html $OUTPUT_DIR/Logs/$DATE"_execution_timeline.html"
