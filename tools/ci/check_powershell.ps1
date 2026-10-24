param(
    [string[]]$Paths,
    [ValidateSet('Auto', 'Windows51', 'PowerShell7')]
    [string]$Compatibility = 'Auto'
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).ProviderPath
$requiredVersion = '1.25.0'
$module = Get-Module -ListAvailable -Name PSScriptAnalyzer |
    Where-Object { $_.Version -eq [version]$requiredVersion } |
    Select-Object -First 1
if ($null -eq $module) {
    throw "PSScriptAnalyzer $requiredVersion is required. Install-Module PSScriptAnalyzer -RequiredVersion $requiredVersion -Scope CurrentUser"
}
Import-Module $module.Path -Force -ErrorAction Stop
$settingsPath = Join-Path $PSScriptRoot 'PSScriptAnalyzerSettings.psd1'
$windowsSettings = Import-PowerShellDataFile -LiteralPath $settingsPath
$coreSettings = Import-PowerShellDataFile -LiteralPath $settingsPath
$coreSettings.Rules.PSUseCompatibleSyntax.TargetVersions = @('7.0')
# The plugin builder runs on PS7 on Linux, Windows and macOS. Desktop-only
# Windows profiles would incorrectly reject its RuntimeInformation APIs.
$coreProfiles = @(
    # This bundled Windows profile avoids the duplicate, truncated OSPlatform
    # type entry present in the 17763 PS7 snapshot; no type rule is suppressed.
    'win-8_x64_10.0.14393.0_7.0.0_x64_3.1.2_core',
    'ubuntu_x64_18.04_7.0.0_x64_3.1.2_core'
)
$coreSettings.Rules.PSUseCompatibleCommands.TargetProfiles = $coreProfiles
$coreSettings.Rules.PSUseCompatibleTypes.TargetProfiles = $coreProfiles
if (!$Paths) {
    $Paths = @('tools', 'tests')
}
$files = @()
foreach ($selected in $Paths) {
    $candidate = if ([IO.Path]::IsPathRooted($selected)) { $selected } else { Join-Path $repoRoot $selected }
    if (!(Test-Path -LiteralPath $candidate)) {
        throw "Selected path does not exist: $selected"
    }
    $candidate = (Resolve-Path -LiteralPath $candidate).ProviderPath
    $prefix = $repoRoot + [IO.Path]::DirectorySeparatorChar
    $pathComparison = if ([IO.Path]::DirectorySeparatorChar -eq '\') {
        [StringComparison]::OrdinalIgnoreCase
    } else {
        [StringComparison]::Ordinal
    }
    if (!$candidate.StartsWith($prefix, $pathComparison)) {
        throw "Selected paths must stay inside this repository: $selected"
    }
    # Resolve-Path normalizes provider paths but does not resolve every directory
    # link target. Reject links anywhere below this checkout before recursion.
    $ancestor = $candidate
    while (![string]::Equals($ancestor, $repoRoot, $pathComparison)) {
        $ancestorItem = Get-Item -LiteralPath $ancestor -Force
        if ($ancestorItem.Attributes -band [IO.FileAttributes]::ReparsePoint) {
            throw "Selected paths must not traverse symbolic links or reparse points: $selected"
        }
        $ancestor = Split-Path -Parent $ancestor
    }
    $relative = $candidate.Substring($prefix.Length).Replace('\', '/')
    if ($relative -notmatch '^(tools|tests)(/|$)') {
        throw "Select maintained paths under tools/ or tests/: $selected"
    }
    $items = if (Test-Path -LiteralPath $candidate -PathType Container) {
        @(Get-ChildItem -LiteralPath $candidate -Recurse -File)
    } else {
        @(Get-Item -LiteralPath $candidate)
    }
    foreach ($item in $items) {
        $relativeFile = $item.FullName.Substring($prefix.Length).Replace('\', '/')
        if ($item.Extension -in @('.ps1', '.psm1') -and
            $relativeFile -notmatch '(^|/)(archives|build|_dist|GH_logs|vendor|generated|\.venv|venv|__pycache__)(/|$)' -and
            !($item.Attributes -band [IO.FileAttributes]::ReparsePoint)) {
            $files += $item.FullName
        }
    }
}
$files = @($files | Sort-Object -Unique)
if ($files.Count -eq 0) {
    throw 'Zero eligible maintained PowerShell files selected (check exclusions).'
}
$diagnostics = @()
foreach ($file in $files) {
    $settings = $windowsSettings
    if ($Compatibility -eq 'PowerShell7' -or
        ($Compatibility -eq 'Auto' -and [IO.Path]::GetFileName($file) -eq 'plugin_build_check.ps1')) {
        $settings = $coreSettings
    }
    $diagnostics += @(Invoke-ScriptAnalyzer -Path $file -Settings $settings -ErrorAction Stop)
}
if ($diagnostics.Count -gt 0) {
    $diagnostics | Format-Table -AutoSize ScriptName, Line, Severity, RuleName, Message | Out-String -Width 240 | Write-Output
    throw "PowerShell analysis failed with $($diagnostics.Count) diagnostic(s); all selected severities fail."
}
Write-Output "PowerShell analysis passed ($($files.Count) files; PSScriptAnalyzer $requiredVersion)."
