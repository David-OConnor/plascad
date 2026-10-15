$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest
$appName = "plascad"

Push-Location -LiteralPath $PSScriptRoot
try {
    # The version is defined in one place only: Cargo.toml's `[package] version`. Read it from there,
    # so a release never ships under a version that disagrees with the one the app reports.
    $manifest = Get-Content -LiteralPath "Cargo.toml" -Raw
    $packageMatch = [regex]::Match($manifest, '(?ms)^\[package\]\s*\r?\n(?<body>.*?)(?=^\[|\z)')
    $versionMatch = [regex]::Match($packageMatch.Groups['body'].Value, '(?m)^\s*version\s*=\s*"([^"]+)"')
    if (-not $versionMatch.Success) { throw "Could not read the version from Cargo.toml" }
    # Archive names use underscores, e.g. 0.1.0 -> 0_1_0.
    $version = $versionMatch.Groups[1].Value -replace '\.', '_'
    $outDir = "target\release"
    $exe    = Join-Path $outDir "${appName}.exe"
    $readme = "README.md"
    $setup = "install_scripts/setup_${appName}.ps1"
    $setupLauncher = "install_scripts/setup_${appName}.bat"

    cargo build --release
    if ($LASTEXITCODE -ne 0) {
        throw "cargo build failed with exit code $LASTEXITCODE."
    }

    $payload = @($exe, $readme, $setup, $setupLauncher)
    foreach ($file in $payload) {
        if (-not (Test-Path -LiteralPath $file -PathType Leaf)) {
            throw "Release file $file was not found."
        }
    }

    $archive = "${appName}_${version}_win.zip"
    if (Test-Path -LiteralPath $archive) { Remove-Item -LiteralPath $archive -Force }
    Compress-Archive -LiteralPath $payload -DestinationPath $archive -Force

    # Build and package the Linux release too, in the Ubuntu 22.04 WSL instance. 22.04 specifically,
    # since it sets the glibc version the Linux binary requires; building on a newer Ubuntu would
    # produce a binary that refuses to start on older distros.
    $wslDistro = "Ubuntu-22.04"
    # `wsl -l` writes UTF-16, which PowerShell reads with embedded nulls; strip them before comparing.
    $wslDistros = (wsl.exe -l -q) | ForEach-Object { ($_ -replace "`0", "").Trim() } | Where-Object { $_ }
    if ($LASTEXITCODE -ne 0) {
        throw "Could not list WSL instances (exit code $LASTEXITCODE)."
    }
    if ($wslDistros -notcontains $wslDistro) {
        throw "The WSL instance '$wslDistro' was not found, so the Linux release was not built. It must be Ubuntu 22.04. Installed instances: $($wslDistros -join ', ')"
    }

    # `--exec` skips the distro's shell, so os-release is read as-is and parsed here.
    $wslOsRelease = wsl.exe -d $wslDistro --exec cat /etc/os-release
    if ($LASTEXITCODE -ne 0) {
        throw "Could not read /etc/os-release in '$wslDistro' (exit code $LASTEXITCODE)."
    }
    $wslVersion = $wslOsRelease | Select-String -Pattern '^VERSION_ID="?([^"]+)"?' |
        ForEach-Object { $_.Matches[0].Groups[1].Value } | Select-Object -First 1
    if ($wslVersion -ne "22.04") {
        throw "The WSL instance '$wslDistro' is not Ubuntu 22.04 (VERSION_ID: '$wslVersion'), so the Linux release was not built."
    }

    # A login shell, so ~/.profile puts cargo on the PATH.
    wsl.exe -d $wslDistro --cd $PSScriptRoot -- bash -l ./make_release.sh
    if ($LASTEXITCODE -ne 0) {
        throw "make_release.sh failed in WSL ($wslDistro) with exit code $LASTEXITCODE."
    }
} finally {
    Pop-Location
}
