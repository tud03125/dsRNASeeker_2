#!/usr/bin/env bash
set -euo pipefail
exec conda run -n sprint_env --no-capture-output python /rs01/home/levinm/dsRNASeeker_2/tools/SPRINT/utilities/getA2I.py "$@"
