@echo off
cd /d "C:\xrt\src"
echo Processing fit files...
powershell -ExecutionPolicy Bypass -NoProfile -Command "& 'C:\Tools\python\venv\callisto_env\Scripts\Activate.ps1'; python -m proc.read_fitfiles_freq"
echo Generating diagnostics plot...
powershell -ExecutionPolicy Bypass -NoProfile -Command "& 'C:\Tools\python\venv\callisto_env\Scripts\Activate.ps1'; python -m proc.plot_solarflux_freq -remove_pkl"
echo All done!