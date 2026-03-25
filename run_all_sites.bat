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
set "OUTDIR=%DIR%\benchmarks"
set "BIN=%DIR%\%BUILD_PATH%"
set "GUESS=%BIN%\guesscmd.exe"

:: Define sites to run
set "SITES_FILE=%OUTDIR%\EC_grass_sites.txt"

if not exist "%SITES_FILE%" (
    echo ERROR: Site list file not found: "%SITES_FILE%"
    echo Please create benchmarks\ozflux\sites.txt with one site per line.
    goto :fail
)

:: Start all simulations
for /f "usebackq delims=" %%S in ("%SITES_FILE%") do (
    set "SITE_DIR=%OUTDIR%\%%S"
    if not exist "!SITE_DIR!" (
        echo Creating site directory "!SITE_DIR!"
        mkdir "!SITE_DIR!"
    ) else (
        pushd "!SITE_DIR!"
        if not exist "%%S.ins" (
            echo WARNING: Missing INS file in "!SITE_DIR!": "%%S.ins"
        ) else (
            start /b "" "%GUESS%" -input nc "%%S.ins"
            rem Optional: log outputs per site
            rem start /b "" "%GUESS%" -input nc "%%S.ins" > "%%S.out.log" 2> "%%S.err.log"
        )
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
