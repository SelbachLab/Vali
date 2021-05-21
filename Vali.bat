set mypath=%cd%
cd %mypath%
set "path1=%mypath%\R-Portable\bin\x64\Rscript.exe start-win.R %mypath%"
@echo %path1%
timeout 3
%path1%
start %path1%

timeout 10
