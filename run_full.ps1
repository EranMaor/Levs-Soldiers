# Full mutual information, then full block error rate.
# BLER starts only after the MI process exits. Each part has its own log.
Set-Location $PSScriptRoot
& "$PSScriptRoot\run_full_mi.ps1"
Write-Host "MI process ended. Starting full BLER."
& "$PSScriptRoot\run_full_bler.ps1"
