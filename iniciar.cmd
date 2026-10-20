@echo off
setlocal
cd /d "%~dp0"
set "HELICE_PYTHON=%USERPROFILE%\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
if exist "%HELICE_PYTHON%" (
  "%HELICE_PYTHON%" server.py %*
) else (
  py -3 server.py %*
)
if errorlevel 1 (
  echo.
  echo Nao foi possivel iniciar. Confira o erro acima e instale Python 3.11 ou superior se necessario.
  pause
)
