@echo off
CHCP 65001
echo [INFO] Starting GlycoLipid-Insight Analysis...
echo [INFO] Target File: data/raw/Ginsenoside.csv

REM Use the specific python interpreter
REM limit 0 means process all records
"C:\Users\Daniel Yang\anaconda3\envs\chem_ai\python.exe" scripts/run_pipeline.py --source local --keyword "data/raw/Ginsenoside.csv" --limit 0

if errorlevel 1 (
    echo [ERROR] Pipeline failed!
    pause
    exit /b %errorlevel%
)

echo [SUCCESS] Pipeline finished. Check reports/unified_data.xlsx and images/Ginsenoside/
echo [INFO] Debugging single sugar visualization...
"C:\Users\Daniel Yang\anaconda3\envs\chem_ai\python.exe" scripts/debug_sugar.py

pause
