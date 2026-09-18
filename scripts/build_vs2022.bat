@echo off
setlocal

set "VSDEVCMD=C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat"
set "CMAKE=C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"

if not exist "%VSDEVCMD%" (
    echo Visual Studio developer command prompt was not found:
    echo %VSDEVCMD%
    exit /b 1
)

if not exist "%CMAKE%" (
    echo Visual Studio bundled CMake was not found:
    echo %CMAKE%
    exit /b 1
)

call "%VSDEVCMD%" -arch=x64
if errorlevel 1 exit /b %errorlevel%

if exist build\CMakeCache.txt (
    "%CMAKE%" -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
    if errorlevel 1 exit /b %errorlevel%
) else (
    "%CMAKE%" -S . -B build -G "Visual Studio 17 2022" -A x64 -DCMAKE_BUILD_TYPE=RelWithDebInfo
    if errorlevel 1 exit /b %errorlevel%
)

"%CMAKE%" --build build --config RelWithDebInfo --parallel
if errorlevel 1 exit /b %errorlevel%

echo.
echo Build completed. Run executables under bin\
