#!/usr/bin/env bash
set -euo pipefail
exec conda run -n sprint_env --no-capture-output sprint "$@"
