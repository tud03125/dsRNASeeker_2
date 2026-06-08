# Installing dsRNASeeker_2 on FCCC HPC

```bash
cd /rs01/home/levinm
unzip dsRNASeeker_2.zip
cd dsRNASeeker_2

# mamba is recommended if available; conda also works but may solve slowly.
mamba env create -f environment_dsRNASeeker_2.yml
# or:
# conda env create -f environment_dsRNASeeker_2.yml

conda activate dsRNASeeker_2
python3 main.py --help
python3 main.py workflow --help
```

If the REDItools executable installed on FCCC is not `REDItoolDnaRna.py`, pass it explicitly with:

```bash
--reditools-exe /path/to/your/REDItoolDnaRna.py
```

If `rmats.py` is not on PATH after activation, pass:

```bash
--rmats-exe /path/to/rmats.py
```

SPRINT is intentionally optional. Add `--run-sprint --sprint-exe /path/to/SPRINT` only after the core workflow is passing.

## REDItools2 note

The environment file does not assume one universal REDItools2 executable name because REDItools2 installations vary. Keep using your existing FCCC REDItools2 installation and point dsRNASeeker_2 to the executable with `--reditools-exe`. If you prefer to skip it during the first workflow test, use `--skip-reditools` or `--precomputed-redit-dir`.
