#!/bin/bash
module load qbic/anaconda2/2.1.0

workflowDir=$(cat wfdir)
#parse using CTDopts and run workflow
python runWorkflow.py $workflowDir
cp wfdir wfdir2