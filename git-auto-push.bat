@echo off
ECHO
echo Git push at %date% %time% >> C:\xrt\git-log.txt
cd /d C:\xrt
"C:\Program Files\Git\cmd\git.exe" add . >> C:\xrt\git-log.txt 2>&1
"C:\Program Files\Git\cmd\git.exe" add .
"C:\Program Files\Git\cmd\git.exe" commit -m "TEST12 Daily automatic update %date%"
"C:\Program Files\Git\cmd\git.exe" push origin master