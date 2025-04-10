@echo off
cd /d "C:\xrt\src"
powershell -ExecutionPolicy Bypass -NoProfile -Command "& 'C:\Tools\python\venv\callisto_env\Scripts\Activate.ps1'; python -m proc.update_daily_flux"