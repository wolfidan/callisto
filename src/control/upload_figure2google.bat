cd /d "C:\xrt\src"
powershell -ExecutionPolicy Bypass -NoProfile -Command "& 'C:\Tools\python\venv\callisto_env\Scripts\Activate.ps1'; python -m control.upload_figure2google"
pause