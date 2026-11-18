# Overnight mutual-information run (paper-scale Fig. 8).
# Do not start this while plots.py is already running in another terminal.
# Output is also written to results\full_mi.log so you can close the window's scrollback.
Set-Location $PSScriptRoot
New-Item -ItemType Directory -Force -Path results | Out-Null
Write-Host "Full MI started $(Get-Date -Format 'yyyy-MM-dd HH:mm'). Log: results\full_mi.log"
python -u plots.py mi --mode full 2>&1 | Tee-Object -FilePath results\full_mi.log
Write-Host "Full MI finished $(Get-Date -Format 'yyyy-MM-dd HH:mm')"
