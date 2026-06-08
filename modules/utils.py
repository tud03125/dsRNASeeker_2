from __future__ import annotations
import os, shutil, subprocess
from pathlib import Path
import pandas as pd


def which_or_none(exe: str) -> str | None:
    return shutil.which(exe)


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def is_nonempty_file(path: str | Path) -> bool:
    p = Path(path)
    return p.exists() and p.is_file() and p.stat().st_size > 0


def run_cmd(
    cmd: list[str],
    env: dict[str, str] | None = None,
    cwd: str | None = None,
    capture: bool = False,
    log_path: str | Path | None = None,
    quiet: bool = False,
) -> subprocess.CompletedProcess:
    """Run a command with optional log redirection.

    quiet=True sends stdout/stderr to log_path when provided, otherwise captures it.
    This keeps the user-facing workflow output at step-level while preserving
    full tool logs for debugging/reproducibility.
    """
    merged = os.environ.copy()
    if env:
        merged.update({k: str(v) for k, v in env.items() if v is not None})

    if log_path is not None:
        lp = Path(log_path)
        lp.parent.mkdir(parents=True, exist_ok=True)
        with lp.open("a", encoding="utf-8") as log:
            log.write("\n# COMMAND\n" + " ".join(map(str, cmd)) + "\n")
            log.flush()
            return subprocess.run(cmd, check=True, env=merged, cwd=cwd, text=True, stdout=log, stderr=subprocess.STDOUT)

    if quiet:
        return subprocess.run(cmd, check=True, env=merged, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    return subprocess.run(cmd, check=True, env=merged, cwd=cwd, text=True, capture_output=capture)


def read_samplesheet(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t')
    required = {'sample_id', 'condition', 'bam_path'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f'Samplesheet missing required columns: {sorted(missing)}')
    return df


def step(message: str) -> None:
    print(f"[dsRNASeeker] {message}", flush=True)
