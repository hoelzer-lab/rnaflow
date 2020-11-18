#!/usr/bin/bash

# copy final reports to output dir

OUTPUT_DIR=$1
REPORTS_DIR=$2

mkdir -p $OUTPUT_DIR/Logs

cp $REPORTS_DIR/execution_report.html $OUTPUT_DIR/Logs
cp $REPORTS_DIR/execution_timeline.html $OUTPUT_DIR/Logs
