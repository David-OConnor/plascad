[CmdletBinding()]
param(
    [string]$SourceDirectory,
    [string]$InstallDirectory
)

# This file installs the application for the current user, and creates a Start menu entry.
#
# Install PlasCAD in %LOCALAPPDATA%\Programs without administrator rights.
# Preferences are saved in that directory, since the Start menu entry starts PlasCAD there.
#
# Run it from the folder you extracted the release into, e.g. by double-clicking
# setup_plascad.bat, or with:
#   powershell -NoProfile -ExecutionPolicy Bypass -File setup_plascad.ps1

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$NAME_UPPER = "PlasCAD"
$NAME = "plascad"
$EXE = "$NAME.exe"
$DESCRIPTION = "Plasmid and primer design"

function Copy-Payload {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Source,
        [Parameter(Mandatory = $true)]
        [string]$Destination
    )

    $exeSource = Join-Path $Source $EXE
    if (-not (Test-Path -LiteralPath $exeSource -PathType Leaf)) {
        throw "$EXE was not found in $Source. Run this script from the folder you extracted the release into."
    }

    if (Get-Process -Name $NAME -ErrorAction SilentlyContinue) {
        throw "$NAME_UPPER is currently running, so its files can't be replaced. Close it, then run this script again."
    }

    New-Item -ItemType Directory -Path $Destination -Force | Out-Null
    if ($exeSource -ne (Join-Path $Destination $EXE)) {
        Copy-Item -LiteralPath $exeSource -Destination $Destination -Force
    }
    Write-Host "Copied $EXE to $Destination."
}

function New-StartMenuEntry {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Destination
    )

    $programs = [Environment]::GetFolderPath("Programs")
    if (-not $programs) {
        throw "The Start menu folder for this user could not be located."
    }
    New-Item -ItemType Directory -Path $programs -Force | Out-Null

    $shortcutPath = Join-Path $programs "$NAME_UPPER.lnk"
    $target = Join-Path $Destination $EXE

    $shell = New-Object -ComObject WScript.Shell
    try {
        $shortcut = $shell.CreateShortcut($shortcutPath)
        $shortcut.TargetPath = $target
        # Start relative paths in the installation directory.
        $shortcut.WorkingDirectory = $Destination
        # The executable embeds its own icon; see build.rs.
        $shortcut.IconLocation = "$target,0"
        $shortcut.Description = $DESCRIPTION
        $shortcut.Save()
    } finally {
        [Runtime.InteropServices.Marshal]::ReleaseComObject($shell) | Out-Null
    }

    Write-Host "Created the Start menu entry at $shortcutPath."
}

# We default this here rather than in the param block: Windows PowerShell 5.1 leaves $PSScriptRoot
# empty while binding parameters when a script is run with -File, as setup_plascad.bat does. It is
# populated by the time the body runs, in every host.
if (-not $SourceDirectory) {
    $SourceDirectory = $PSScriptRoot
}
if (-not $SourceDirectory) {
    $SourceDirectory = (Get-Location).Path
}
if (-not (Test-Path -LiteralPath $SourceDirectory -PathType Container)) {
    throw "The source directory $SourceDirectory does not exist."
}
$SourceDirectory = (Resolve-Path -LiteralPath $SourceDirectory).ProviderPath

if (-not $InstallDirectory) {
    if (-not $env:LOCALAPPDATA) {
        throw "LOCALAPPDATA is not set; provide InstallDirectory explicitly."
    }
    $InstallDirectory = Join-Path $env:LOCALAPPDATA "Programs\$NAME_UPPER"
}
# Shortcuts require an absolute target even when InstallDirectory was supplied as a relative path.
$InstallDirectory = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($InstallDirectory)

Copy-Payload -Source $SourceDirectory -Destination $InstallDirectory
New-StartMenuEntry -Destination $InstallDirectory

Write-Host ""
Write-Host "$NAME_UPPER is installed in $InstallDirectory."
Write-Host "You can launch it from the Start menu (e.g., search `"$NAME`"), and/or pin it to the taskbar."
Write-Host ""
