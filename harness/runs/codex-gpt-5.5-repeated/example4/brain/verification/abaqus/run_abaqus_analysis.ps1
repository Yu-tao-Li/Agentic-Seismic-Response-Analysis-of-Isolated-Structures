$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ScriptDir

python generate_abaqus_input.py *> abaqus_generate.log
$jobFiles = @(
  "isolated_shear_building.com", "isolated_shear_building.dat", "isolated_shear_building.env",
  "isolated_shear_building.ipm", "isolated_shear_building.log", "isolated_shear_building.msg",
  "isolated_shear_building.odb", "isolated_shear_building.odb_f", "isolated_shear_building.prt",
  "isolated_shear_building.sim", "isolated_shear_building.sta"
)
foreach ($file in $jobFiles) {
  if (Test-Path -LiteralPath $file) {
    Remove-Item -LiteralPath $file -Force
  }
}
$ErrorActionPreference = "Continue"
abaqus job=isolated_shear_building input=isolated_shear_building.inp interactive *> abaqus_run.log
$abaqusExit = $LASTEXITCODE
$ErrorActionPreference = "Stop"
if (-not (Test-Path -LiteralPath "isolated_shear_building.sta")) {
  throw "ABAQUS did not create isolated_shear_building.sta; native exit code $abaqusExit"
}
$staText = Get-Content -LiteralPath "isolated_shear_building.sta" -Raw
if ($staText -notmatch "THE ANALYSIS HAS COMPLETED SUCCESSFULLY") {
  throw "ABAQUS analysis did not complete successfully; native exit code $abaqusExit"
}
abaqus python postprocess_abaqus_odb.py *> abaqus_postprocess_stdout.log
