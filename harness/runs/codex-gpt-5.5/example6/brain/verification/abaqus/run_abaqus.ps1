$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $scriptDir

python generate_inp.py

$vsDir = "<visual_studio_buildtools>"
$setvars = "<intel_oneapi>\setvars.bat"
$cmd = @"
set "VS2022INSTALLDIR=$vsDir" && call "$setvars" intel64 vs2022 --force && abaqus job=example6_model input=example6_model.inp user=boucwen_uel.for interactive && abaqus python extract_results.py
"@

cmd.exe /d /s /c $cmd
