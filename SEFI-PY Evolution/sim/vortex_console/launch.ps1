$ErrorActionPreference = 'Stop'
$bundledPython = Join-Path $env:USERPROFILE '.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
if (Test-Path -LiteralPath $bundledPython) { $consolePython = $bundledPython } else { $consolePython = (Get-Command python -ErrorAction Stop).Source }
Set-Location -LiteralPath $PSScriptRoot
Write-Host 'Open http://127.0.0.1:8765 in your browser. Keep this console running; Ctrl+C stops it.'
& $consolePython server.py
