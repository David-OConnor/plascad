#!/usr/bin/env bash
set -euo pipefail

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
app_name="plascad"

# The version is defined in one place only: Cargo.toml's `[package] version`. Read it from there,
# so a release never ships under a version that disagrees with the one the app reports.
# Archive names use underscores, e.g. 0.1.0 -> 0_1_0.
version=$(awk '
    /^\[package\][[:space:]]*$/ { in_package = 1; next }
    /^[[:space:]]*\[/ { in_package = 0 }
    in_package && /^[[:space:]]*version[[:space:]]*=/ {
        split($0, parts, "\""); print parts[2]; exit
    }
' Cargo.toml | tr '.' '_')
if [[ -z "$version" ]]; then
    echo "Could not read the version from Cargo.toml" >&2
    exit 1
fi
outdir=target/release
exe="$outdir/$app_name"
readme="README.md"
setup="install_scripts/setup_${app_name}.sh"
icon="src/resources/icon.png"

cargo build --release

for file in "$exe" "$readme" "$setup" "$icon"; do
    if [[ ! -f "$file" ]]; then
        printf 'Error: release file %s was not found.\n' "$file" >&2
        exit 1
    fi
done
chmod +x "$exe" "$setup"

# The zip takes the setup script as it is on disk, not as git stores it. On Windows, that copy can
# have CRLF endings (core.autocrlf, or an editor), which makes bash fail on the first line for
# whoever installs the release. Strip them here so the archive is right regardless.
sed -i 's/\r$//' "$setup"

# Recreate the archive so files from an earlier release cannot remain in it.
archive="${app_name}_${version}_linux.zip"
rm -f -- "$archive"
zip -j "$archive" "$exe" "$readme" "$setup" "$icon"
