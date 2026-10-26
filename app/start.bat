@echo off
chcp 65001 >nul
cd /d "%~dp0"
python -c "import numpy, scipy" >nul 2>nul
if errorlevel 1 (
  echo 需要 Python 3.10 或更高版本，以及 numpy、scipy。
  echo 安装依赖：python -m pip install -r requirements.txt
  pause
  exit /b 1
)
python -X utf8 server.py
pause
