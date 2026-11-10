@echo off
setlocal
cd /d "%~dp0.."

if not exist ".venv\Scripts\python.exe" (
  echo Project Python environment was not found.
  echo Run the setup instructions in README.md first.
  pause
  exit /b 1
)

".venv\Scripts\python.exe" "ui\build_ui_data.py"
if errorlevel 1 (
  echo Could not refresh the peptide lookup data.
  pause
  exit /b 1
)

start "" "%CD%\ui\index.html"
endlocal
