# Stop the most recently started dpenum process if free physical memory drops below a threshold.
param([int]$MinFreeMB = 400, [string]$Log = "results/logs/watchdog.log")
while ($true) {
  $os = Get-CimInstance Win32_OperatingSystem
  $free = [math]::Round($os.FreePhysicalMemory / 1KB)
  $procs = Get-Process -Name dpenum* -ErrorAction SilentlyContinue | Sort-Object StartTime
  if (-not $procs) { Start-Sleep -Seconds 30; continue }
  if ($free -lt $MinFreeMB -and $procs.Count -gt 1) {
    $victim = $procs[-1]
    "$(Get-Date -Format s) free=${free}MB stopping pid $($victim.Id) (started $($victim.StartTime))" | Out-File -Append -Encoding utf8 $Log
    Stop-Process -Id $victim.Id -Force
  }
  Start-Sleep -Seconds 20
}
