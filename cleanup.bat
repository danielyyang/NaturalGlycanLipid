@echo off
echo =========================================================
echo GlycoNP Project Cleanup Utility
echo ---------------------------------------------------------
echo This script will permanently delete outdated prototype
echo scripts and V11/V12/V13 artifacts to clean the repo.
echo Ensure you have committed important files to git!
echo =========================================================
pause

REM Change to script directory
cd /d "%~dp0"

echo Deleting root level temporary files...
del /q "tmp_check_hex.py" >nul 2>&1
del /q "rematch_log.txt" >nul 2>&1
del /q "start_project.bat" >nul 2>&1

echo Deleting V11/V12/V13 legacy pipelines...
del /q "scripts\run_v11_final.py" >nul 2>&1
del /q "scripts\run_v12_full_pipeline.py" >nul 2>&1
del /q "scripts\run_v13_full_pipeline.py" >nul 2>&1
del /q "scripts\main_pipeline.py" >nul 2>&1
del /q "scripts\run_glyconp_pipeline.py" >nul 2>&1

echo Deleting obsolete chart splitting scripts...
del /q "scripts\generate_saponin_charts_part2.py" >nul 2>&1
del /q "scripts\generate_saponin_charts_part3.py" >nul 2>&1
del /q "scripts\generate_executive_plots.py" >nul 2>&1
del /q "scripts\generate_glyco_landscape.py" >nul 2>&1
del /q "scripts\plot_family_network.py" >nul 2>&1
del /q "scripts\plot_scaffold_taxonomy_map.py" >nul 2>&1
del /q "scripts\plot_umap_chemical_space.py" >nul 2>&1

echo Deleting hallucinated NLP guesser...
del /q "scripts\rescue_generic_sugars_phase2.py" >nul 2>&1

echo Deleting PoC engine implementations...
del /q "scripts\cip_exo_identification_engine.py" >nul 2>&1
del /q "scripts\poc_cip_fingerprint_engine.py" >nul 2>&1
del /q "scripts\poc_virtual_demodify.py" >nul 2>&1

echo Deleting diagnostic and debug scripts...
del /q "scripts\_check_scaffolds.py" >nul 2>&1
del /q "scripts\_debug_scaffolds.py" >nul 2>&1
del /q "scripts\_debug_scaffolds_deep.py" >nul 2>&1
del /q "scripts\compare_v13_v_new.py" >nul 2>&1
del /q "scripts\preview_saponin_columns.py" >nul 2>&1
del /q "scripts\diag_last4.py" >nul 2>&1
del /q "scripts\diag_vitexin.py" >nul 2>&1
del /q "scripts\_check_hex_prop.py" >nul 2>&1
del /q "scripts\_check_saponin_stats.py" >nul 2>&1
del /q "scripts\_check_scaffold_dist.py" >nul 2>&1

echo Deleting outdated Benchmark tier run-scripts...
del /q "scripts\run_benchmark_200_tier_a.py" >nul 2>&1
del /q "scripts\run_benchmark_report.py" >nul 2>&1
del /q "scripts\build_benchmark_150.py" >nul 2>&1
del /q "scripts\build_benchmark_200.py" >nul 2>&1

echo Cleanup complete.
pause
