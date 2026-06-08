#!/usr/bin/env bash
set -euo pipefail

export PYTHONNOUSERSITE=1

exec /rs01/projects/jadezhoulab/tud03125/anaconda3/envs/dott_pipeline_Temple/dsRNASeeker_2/bin/python \
  /rs01/home/levinm/dsRNASeeker_2/tools/REDItools2/src/cineca/reditools.py "$@"