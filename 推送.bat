@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"

where git >nul 2>nul
if errorlevel 1 (
    echo [ERROR] git not found. Install Git for Windows first: https://git-scm.com/download/win
    pause
    exit /b 1
)

if not exist ".git" (
    echo [1/5] git init
    git init
    git branch -M main
)

echo [2/5] git add + commit
git add -A
git -c user.name="kunstk" -c user.email="kunstk@users.noreply.github.com" commit -m "update: His-MS deconvolution" 2>nul

echo [3/5] remote origin
git remote get-url origin >nul 2>nul
if errorlevel 1 (
    git remote add origin https://github.com/dutaodang-crypto/INK.git
)

echo [4/5] git push
git push -u origin main

if errorlevel 1 (
    echo.
    echo [retry] the remote may already have commits, trying rebase...
    git pull --rebase origin main
    git push -u origin main
)

echo.
echo Finished. If a password is requested, use a Personal Access Token:
echo   GitHub - Settings - Developer settings - Personal access tokens (classic) - scope: repo
echo.
pause
