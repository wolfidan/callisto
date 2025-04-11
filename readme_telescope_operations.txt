===================================================
Telescope installation and operation
===================================================

file:///C:/CALLISTO-01/Application/eCallistoManual.pdf

Installation and alignment
(0.1) Install telescope and connect cables (without current on the system) 
(0.2) Switch DiSEqC controller & CALLISTO spectrometer
(0.3) Align telescope: see separate section below

Automatic scanning routine
(1.1) Check that "CaliProc = True" in following configuration file:
      C:\xrt\src\PythonScripts\TrackingSun\configsun.cfg
(1.2) Ev. edit frequency file, e.g. e.g. C:\CALLISTO-01\Application\freq99995.cfg
(1.3) Edit the configuration file "C:\CALLISTO-01\Application\scheduler.cfg":  
(1.4) Open system scheduler and activate "Suntracker" cron job
(1.5) Open the application software ‘callisto.exe’: set automatic, 
      ev. open lightcurve (shows data as soon as the measurements start)

Remarks:
- Measurements are automatically started/stopped according to "C:\CALLISTO-01\Application\scheduler.cfg".
- callisto.exe continuously measures, results are stored in "C:\xrt\output\data\raw\FITfiles" (one every 15 minutes)
- at the same time, system scheduler calls script "C:\xrt\src\PythonScripts\TrackingSun\sunpos_AZI_ELE.py":
  move antenna to sun, or ground/sky (provided that "CaliProc = True" in configsun.cfg)

----------------------------------------
Telescope alignment: manual sun tracking
----------------------------------------

Manual interaction with satellite rotors: file:///C:/xrt/docs/TelescopSetup-V20230914.pdf (in particular, p. 9)
(0.1) open putty, load and open DiSEqC session (9600 baud, 8 data bits, 1 stop bit, no parity and flow control None)
(0.2) in the putty shell, type: -h
(0.3) now we can manually move the antenna with commands like: azi6, ele-1, ele0, azi-30, azi0, ...
      Here we want to use azi0, ele0. 
      Then you may mechanically fine adjust elevation using a balance.
(0.4) Close putty (when Putty is connected to rotor controllers, then other scripts have no access)
(0.5) Deactivate system scheduler application "suntracker"

(1) modify C:\xrt\src\PythonScripts\TrackingSun\configsun.cfg (CaliProc = False)
(2) Open application software "callisto.exe": set manual, select frequency file (e.g. C:\CALLISTO-01\Application\freq99995.cfg), 
    plot lightcurve
(3) Open PowerShell do the following:
    C:\Tools\python\venv\callisto_env\Scripts\Activate.ps1
    cd C:\xrt\src\PythonScripts\TrackingSun
    python .\sunpos_AZI_ELE.py
    --> This points the telescope to the computed position of the sun,
        but probably the alignment is not yet perfect.
    --> Edit the offset values aziref & eleref in following configuration file:
        C:\xrt\src\PythonScripts\TrackingSun\configsun.cfg
    --> re-run the script and test whether the alignment improved 
        (lightcurve value increase, bright spot in the middle of the LNB)
    Once the optimal reference values are found, we leave the environment with the following command:
    deactivate
(4) modify C:\xrt\src\PythonScripts\TrackingSun\configsun.cfg (CaliProc = True)
(5) activate system scheduler application "suntracker"

-------------------------------------------------
Processing the raw measurements --> solarflux
-------------------------------------------------

1. Open the PowerShell.

2. Activate the virtual environment:

   C:\Tools\python\venv\callisto_env\Scripts\Activate.ps1

3. Change to correct folder:

   cd C:\xrt\src

4. Create solarflux values from raw data.

   python -m proc.read_fitfiles
   python -m proc.read_fitfiles -verbose (for printing info to screen)

   Remarks: 
   - By default, the above script reads all the FIT-files "C:\xrt\output\data\raw\FITfiles".
   - For each day, the processed data is stored in a dedicated folder.
     For day 20240510 the file is the following: "C:\xrt\output\data\proc\20240510\solarflux.csv".
   - The computation of the values is not yet fully correct + relies on definitions in the script:
     ambient temperature, frequency range, ... Use e.g. VScode to edit the code.
   - Among all the measurements, Christian suggested to use only the highest 3%. 
     Their average value is used to compute the solarflux estimate: "Ssun_minus_Sback_sfu".

5. Visualize the processed data (e.g. for day 20240510)
    
    cd C:\xrt\src\
    python -m proc.plot_solarflux -day 20240510
    python -m proc.plot_solarflux -day 20240510 -store (for storing the figure)
   
5. Deactivate the virtual environment:
   
   deactivate

-------------------------------------------------
Remarks about organization of scripts & data
-------------------------------------------------

"C:\CALLISTO-01"
By default (I have not found where to change that), the e-CALLISTO software is installed to "C:\CALLISTO-01"

"C:\Tools"
By my choice, additional software is installed for the use of e-CALLISTO is installed in "C:\Tools" (e.g. Arduino, VS Code, ...)

"C:\xrt\"
I created that folder for not saving everything into the above-introduced installation directory.
- C:\xrt\output\log\LogFiles
  Contains the logging files of Christian's software.
- C:\xrt\output\data\raw\FITfiles
  The measurements from the spectrometer are stored here.
- C:\xrt\output\data\proc
  Here we store any data that we generate on the basis of the raw data.

-------------------------------------------------
Additional comments
-------------------------------------------------

- If you don’t want to work with the stored scheduler, 
  then in the GUI press the radio-button labeled ‘manual’.
- The system may run with its internal clock only. But we recommend to automatically 
  synchronizing PC-clock via network to a standard UT-atomic-timing system:
  e.g. NTP Time Server Monitor by Meinberg (I set ut up to launch automatically at start-up)
  --> I set up the system and it should automatically launch at computer start-up.
- Code execution with Python. I set up a virtual environment, within which the code must be executed.
  We eo not install the modules or packages directly to your "system install" of Python: Modules and packages 
  sometimes conflict with each other and with the version of Python you have installed on your system. 
  If there is a compatibility problem, it can cause instability or bugs when you try to use Python. 
  You should use Python's built-in virtual environments instead:
  Each virtual environment can have its own Python version, separate packages and modules, and other variables. 
  That lets you keep the dependencies for each project separate from each other and from your system installation. 
  This ensures that compatibility problems won't affect the primary Python installation on your PC, 
  and that it doesn't become a bloated mess of extra packages and modules.

Possible issues:
- access denied error (e.g. already happened when putty was connected before running the scripts). Resolved when restarting the computer
  https://windowsreport.com/unable-open-serial-port/#:~:text=How%20do%20I%20fix%20unable%20to%20open%20serial,that%20the%20system%20is%20up%20to%20date%20
- check COM ports on windows: Device manager

=====================================================================
Software installation for CALLISTO spectrometer and DiSEqC controller
=====================================================================

-------------------------------------------------
Install software related to CALLISTO spectrometer
- https://e-callisto.org/Software/CallistoInstallerGuide.pdf
- https://www.e-callisto.org/Callisto_DataStatus_Vwww.pdf
- https://e-callisto.org/Software/Callisto-Software.html
-------------------------------------------------

https://www.java.com/de/download/manual.jsp

ssfree (windows, vs. linux: crontab)
https://www.splinterware.com/download/index.html

https://anydesk.com/en/downloads/thank-you?dv=win_exe

https://www.activestate.com/products/perl/
https://platform.activestate.com/xrt-telescope/Perl-5.36.3-Windows/distributions

https://www.reeve.com/Documents/Articles%20Papers/Reeve_NTP-MeinMon_Install.pdf
Create folder C:\Tools\
https://www.tenforums.com/tutorials/135522-add-remove-internet-time-servers-windows.html
https://www.meinbergglobal.com/english/sw/ntp.htm#ntp_stable
https://www.metas.ch/metas/it/home/fabe/zeit-und-frequenz/time-dissemination.html
For deleting file as admin: Start -> Windows PowerShell -> run as administrator -> net user ntp /delete

https://www.meinbergglobal.com/english/sw/ntp-server-monitor.htm
I cannot find the option "Run this program as administrator" here: 
rightclick icon -> properties -> compatibility mode

-------------------------------------------------
Install software related to DiSEqC controller
- https://e-callisto.org/Hardware/Diseqc/SARA%20Diseq_V01.pdf
-------------------------------------------------

https://www.arduino.cc/en/software/download-thank-you
https://www.arduino.cc/en/Guide/DriverInstallation

https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html
https://phoenixnap.com/kb/install-putty-on-windows#:~:text=1%20Installing%20Putty%20on%20Windows%202%20Step%201%3A,official%20website%3Ahttps%3A%2F%2Fwww.chiark.greenend.%203%20Step%202%3A%20Configuration%20and%20Installation

---------------------------------------
Prepare PUTTY session for rotor controller
---------------------------------------

Putty --> Session --> Serial --> COM4
Terminal --> Line discipline options --> Local echo: "Force on" (display text I write)
Terminal --> Line discipline options --> Local line editing: "Force on" (only send text upon pressing Enter)

Select the session DiSEqC I created. Write -h in the terminal and click enter.
- ele (+ upwards, - downwards):
  ele69 --> black sign is about at 63° (higher numbers have no effect)
  ele0 --> black sign at about -7°
  ele69 --> black sign at about -78° (lowe numbers have no effect)
- azi (+ counterclockwise, - clockwise): [-69°, 69°]
- azi0, ele0 when coming from negative numbers: antenna first goes to some positive number and then drives back down to zero

---------------------------------------
Additional software
---------------------------------------

- Install VSCode: https://code.visualstudio.com/docs/?dv=win64user

-------------------------------------------------
Install python on windows
-------------------------------------------------

Install virtual environment
- https://www.howtogeek.com/197947/how-to-install-python-on-windows/
  https://www.python.org/downloads/windows/
  https://docs.python.org/3/library/venv.html
- Do not install the modules or packages directly to your "system install" of Python. Modules and packages 
  sometimes conflict with each other and with the version of Python you have installed on your system. 
  If there is a compatibility problem, it can cause instability or bugs when you try to use Python. 
  You should use Python's built-in virtual environments instead.
- Python lets you create a small virtual environment for each project that you're working on. 
  Each virtual environment can have its own Python version, separate packages and modules, and other variables. 
  That lets you keep the dependencies for each project separate from each other and from your system installation. 
  This ensures that compatibility problems won't affect the primary Python installation on your PC, 
  and that it doesn't become a bloated mess of extra packages and modules.
- venv, or integrated development environment (IDE like PyCharm, ...). 
  IDEs provide all sorts of helpful features if you're coding, and in the case of Python, 
  most include tools to create and manage virtual environments for your projects.
- When a Python interpreter is running from a virtual environment, sys.prefix and sys.exec_prefix point to the 
  directories of the virtual environment, whereas sys.base_prefix and sys.base_exec_prefix point to those of the 
  base Python used to create the environment. It is sufficient to check sys.prefix != sys.base_prefix to determine 
  if the current interpreter is running from a virtual environment:
  --> python -m venv C:\Tools\python\venv\callisto_env

Activate the virtual environment
- cmd.exe vs PowerShell: https://www.varonis.com/blog/powershell-vs-cmd
  In cmd.exe use the following:
  C:\Tools\python\venv\callisto_env\Scripts\activate.bat
  In PowerShell the following:
  C:\Tools\python\venv\callisto_env\Scripts\Activate.ps1

Installations
- install pyserial instead of serial (imp no longer supported in python 3.12)
  https://pyserial.readthedocs.io/en/latest/pyserial_api.html

where python (does not find the version installed in virtual environment)
- https://stackoverflow.com/questions/18713086/virtualenv-wont-activate-on-windows
  https://learn.microsoft.com/en-us/powershell/module/microsoft.powershell.core/about/about_execution_policies?view=powershell-7.4
  Get-ExecutionPolicy -List
  Set-ExecutionPolicy -ExecutionPolicy RemoteSigned

On the e-CALLISTO webpage, there are template python scripts to read the FITfiles.
For now, I copied some of the files to C:\CALLISTO-01, but I think it is best to 
revised code+structure for our purposes. Create folder with new code:
- C:\xrt\src\
- Instead of using Anaconda + integrated development environment, 
  create virtual environment environment and use command prompt.