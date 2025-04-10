@echo off & :: deactivate default printing of comments
rem comments are made with "rem", inline comments with "& ::"
rem execute sun tracker according to location, date and time (from venv)
rem Philipp Schmid, 2024-05-08

set original_dir=%CD%
set venv_root_dir="C:\Tools\python\venv\callisto_env"
set src_dir="C:\xrt\src\PythonScripts\TrackingSun"
call %venv_root_dir%\Scripts\activate.bat
cd %src_dir% & :: script calls configfile in same directory
python sunpos_AZI_ELE.py
cd %original_dir%
call %venv_root_dir%\Scripts\deactivate.bat
rem pause
timeout 1 & :: needed to wait for antenna positioning?
exit /B 1 & :: needed??