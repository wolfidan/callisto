@echo off
echo Starting git operations at %date% %time% >> C:\xrt\git-log.txt
cd /d C:\xrt
"C:\Program Files\Git\cmd\git.exe" add . >> C:\xrt\git-log.txt 2>&1
"C:\Program Files\Git\cmd\git.exe" commit -m "TEST3 Daily automatic update %date%" >> C:\xrt\git-log.txt 2>&1
"C:\Program Files\Git\cmd\git.exe" push origin master >> C:\xrt\git-log.txt 2>&1
echo Finished at %date% %time% >> C:\xrt\git-log.txt