# Run a program and report its peak working set (sampled every 100 ms) and
# wall time. Usage: powershell -File peakmem.ps1 <exe> <args...>
# (plain $args, not a param block, so engine flags like -o and -n pass through)
$Exe = (Resolve-Path $args[0]).Path
$Rest = $args[1..($args.Count - 1)]
$psi = New-Object System.Diagnostics.ProcessStartInfo $Exe
$psi.Arguments = ($Rest | ForEach-Object { if ($_ -match '\s') { '"' + $_ + '"' } else { $_ } }) -join ' '
$psi.UseShellExecute = $false
$psi.RedirectStandardError = $true
$sw = [Diagnostics.Stopwatch]::StartNew()
$psi.WorkingDirectory = (Get-Location).Path
$p = [Diagnostics.Process]::Start($psi)
if (-not $p) { exit 1 }
$errTask = $p.StandardError.ReadToEndAsync()
$peak = 0
while (-not $p.HasExited) {
    try { $p.Refresh(); if ($p.PeakWorkingSet64 -gt $peak) { $peak = $p.PeakWorkingSet64 } } catch {}
    Start-Sleep -Milliseconds 100
}
$p.WaitForExit()
$sw.Stop()
$err = $errTask.Result
($err -split "`n" | Where-Object { $_ -match '^done|^prime' }) -join "`n"
"peak working set MiB {0:N1}  wall {1:N1}s  exit {2}" -f ($peak / 1MB), $sw.Elapsed.TotalSeconds, $p.ExitCode
