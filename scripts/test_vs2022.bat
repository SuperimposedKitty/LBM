@echo off
setlocal

set "CTEST=C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\ctest.exe"

if not exist "%CTEST%" (
    echo Visual Studio bundled CTest was not found:
    echo %CTEST%
    exit /b 1
)

if not exist build\CTestTestfile.cmake (
    echo Test configuration was not found. Run scripts\build_vs2022.bat first.
    exit /b 1
)

"%CTEST%" --test-dir build -C RelWithDebInfo --output-on-failure
exit /b %errorlevel%
