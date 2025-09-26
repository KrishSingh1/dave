@echo off
setlocal enabledelayedexpansion

:: Get the directory containing this script
set "DIR=%~dp0"
cd "%DIR%"

:: Build the model (similar to the bash script)
set "BUILD_PATH=build"
cmake ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DCMAKE_CXX_STANDARD=14 ^
    -DCMAKE_CXX_STANDARD_REQUIRED=ON ^
    -DCMAKE_CXX_EXTENSIONS=ON ^
    -B ^
    "%BUILD_PATH%"

cmake --build "%BUILD_PATH%"

:: Get paths
set "OZFLUX=%DIR%\benchmarks\ozflux"
set "BIN=%DIR%\%BUILD_PATH%"
set "GUESS=%BIN%\guesscmd.exe"

:: Define sites (same list as in bash script)
set "SITES=DryRiver Samford SturtPlains"

:: Start all simulations
for %%S in (%SITES%) do (
    set "SITE_DIR=%OZFLUX%\%%S"
    if not exist "!SITE_DIR!" (
        echo WARNING: Site '%%S' does not appear to exist. Expected: '!SITE_DIR!'
    ) else (
        pushd "!SITE_DIR!"
        start /b "" "%GUESS%" -input nc "%%S.ins"
        popd
    )
)

:: Wait for all background processes to complete
echo Waiting for all simulations to complete...
:wait_loop
timeout /t 5 /nobreak > nul
tasklist | find /i "guesscmd.exe" > nul
if not errorlevel 1 goto wait_loop

echo All simulations completed.

pause