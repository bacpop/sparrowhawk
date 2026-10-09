#!/usr/bin/env python3
"""Build a local wasm-bindgen web package, without invoking the website build.

Python 3.9+ and nightly Rust with rust-src are required. Example:
    python3 scripts/build_wasm64.py --crate rust/ska.rust --module ska
"""

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WASM64_MAX_MEMORY_BYTES = 16 * 1024 * 1024 * 1024



def run(command: list[str]) -> None:
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)



def build(crate: Path, module: str) -> None:
    crate = (ROOT / crate).resolve()
    if not crate.is_relative_to(ROOT):
        raise ValueError("The crate must be inside the sparrowhawk folder.")

    manifest = crate / "Cargo.toml"
    # Resolving the dependencies here (rather than reading Cargo.lock) lets the build work without a
    # committed lockfile, as in CI; cargo writes a local one if it is missing.
    metadata = json.loads(subprocess.check_output([
        "cargo", "metadata", "--format-version", "1",
        "--filter-platform", "wasm64-unknown-unknown",
        "--manifest-path", str(manifest),
    ], cwd=ROOT, text=True))

    package = next(p for p in metadata["packages"] if Path(p["manifest_path"]).resolve() == manifest)
    libraries = [t for t in package["targets"] if "cdylib" in t["crate_types"]]

    if len(libraries) != 1:
        raise ValueError("Expected one cdylib target in the crate manifest.")

    resolved = {node["id"] for node in metadata["resolve"]["nodes"]}
    versions = {p["version"] for p in metadata["packages"]
                if p["name"] == "wasm-bindgen" and p["id"] in resolved}

    if len(versions) != 1:
        raise ValueError("Expected exactly one resolved wasm-bindgen version.")
    
    # Checking versions of wasm-bindgen
    version = versions.pop()
    if tuple(int(part) for part in version.split(".")[:3]) < (0, 2, 120):
        raise ValueError(f"wasm-bindgen {version} predates wasm64 support; update the crate first.")
    
    components = subprocess.check_output([
        "rustup", "component", "list", "--toolchain", "nightly", "--installed",
    ], cwd=ROOT, text=True)
    
    if "rust-src" not in components.splitlines():
        raise ValueError("Install the prerequisite: rustup component add rust-src --toolchain nightly")

    target = ROOT / ".wasm64-build" / module

    # The actual build 
    run([
        "cargo", "+nightly", "build", "--lib", "--release",
        "--manifest-path", str(manifest), "--target", "wasm64-unknown-unknown",
        "-Z", "build-std=std,panic_abort", "--target-dir", str(target),
        "--config",
        f'target.wasm64-unknown-unknown.rustflags=["-C", "link-arg=--max-memory={WASM64_MAX_MEMORY_BYTES}"]',
    ])

    # Now, the binding. THis will check (for local runs) and re-install wasm-bindgen-cli if needed.
    wasm = target / "wasm64-unknown-unknown" / "release" / (libraries[0]["name"] + ".wasm")
    tools = ROOT / ".wasm64-tools" / version
    bindgen = tools / "bin" / "wasm-bindgen"
    if not bindgen.exists():
        run([
            "cargo", "install", "wasm-bindgen-cli", "--locked", "--version", version,
            "--root", str(tools), "--target-dir", str(tools / "target"),
        ])

    actual = subprocess.check_output([str(bindgen), "--version"], text=True).strip()
    
    if actual != f"wasm-bindgen {version}":
        raise ValueError(f"Unexpected CLI version: {actual}; expected {version}.")
    
    staging = target / "web-package"
    staging.mkdir(parents=True, exist_ok=True)
    
    run([str(bindgen), "--target", "web", "--out-name", "index", "--out-dir", str(staging), str(wasm)])
    
    # Checking if we got it!
    if not all((staging / name).is_file() for name in ["index.js", "index_bg.wasm"]):
        raise ValueError("Binding generation did not produce the expected web package.")
    
    # Finally, copying 
    output = ROOT / "www" / "public" / "pkg_wasm64" / module
    shutil.copytree(staging, output, dirs_exist_ok=True)
    print(f"Built {module}: {output}", flush=True)



def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crate", type=Path, required=True)
    parser.add_argument("--module", required=True)
    args = parser.parse_args()
    try:
        build(args.crate, args.module)
    except (ValueError, OSError, subprocess.CalledProcessError, StopIteration) as error:
        print(f"wasm64 build failed: {error}", file=sys.stderr)
        return 1
    return 0



if __name__ == "__main__":
    sys.exit(main())
