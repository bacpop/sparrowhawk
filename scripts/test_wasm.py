#!/usr/bin/env python3
"""Build each crate's browser packages and run its WASM checks under Node: `npm run test:wasm`.

Every check is listed in TESTS. Its wasm32 package is built as a web target into .wasm-test/<module>/wasm32
(the page's own wasm32 builds are bundler targets, which Node cannot load), and a wasm64 package, where the
site has one, by build_wasm64.py into www/public/pkg_wasm64/<module>/. Arguments after `--` go to every
check, e.g. `npm run test:wasm -- --only sparrowhawk -- --verbose`.
"""

import argparse
import subprocess
import sys
from dataclasses import dataclass

from build_wasm64 import ROOT, build as build_wasm64


@dataclass(frozen=True)
class Test:
    """A crate (relative to ROOT), its check (relative to ROOT), the wasm-pack options vue.config.js gives
    it, whether the site builds it for wasm64, and the other modules whose packages its check also needs."""
    crate: str
    check: str
    wasm_pack_args: tuple = ()
    wasm64: bool = False
    uses: tuple = ()


# Keyed by module, as build_wasm64.py names it.
TESTS = {
    "sparrowhawk": Test("rust/sparrowhawk-asm", "tests/wasm/assembler.mjs", wasm64=True),
    "ska": Test("rust/ska.rust", "tests/wasm/ska.mjs", wasm64=True),
    "sketchlib": Test("rust/sketchlib.rust", "tests/wasm/sketchlib.mjs"),
    "amr": Test("rust/sparrowhawk-amr", "tests/wasm/amr.mjs"),
    "orphos": Test("rust/orphos-bridge", "tests/wasm/genecalling.mjs", uses=("amr",)),
    "esm": Test("rust/esm-bridge", "tests/wasm/esm.mjs"),
    "mandrake": Test("rust/mandrake", "tests/wasm/mandrake.mjs",
                     ("--no-default-features", "--features", "wasm-sketchlib")),
}
WASM_TEST = ROOT / ".wasm-test"


def wasm32_package(module: str):
    return WASM_TEST / module / "wasm32"


def build_wasm32(module: str) -> None:
    test = TESTS[module]
    # From the crate's folder, so its .cargo/config.toml (SIMD) applies, as in the page's build.
    subprocess.run(["wasm-pack", "build", "--target", "web", "--release", "--out-dir", str(wasm32_package(module)),
                    "--out-name", "index", *test.wasm_pack_args], cwd=ROOT / test.crate, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", choices=sorted(TESTS), action="append", help="run only this check")
    parser.add_argument("--no-build", action="store_true", help="reuse the packages of the last run")
    parser.add_argument("--update-golden", action="store_true", help="record new golden values")
    parser.add_argument("check_args", nargs=argparse.REMAINDER, help="after --: passed to every check")
    args = parser.parse_args()
    modules = args.only or list(TESTS)
    if not args.no_build:
        for module in dict.fromkeys(m for name in modules for m in (name, *TESTS[name].uses)):
            build_wasm32(module)
            if TESTS[module].wasm64:
                build_wasm64(TESTS[module].crate, module)
    failed = []
    for module in modules:
        test = TESTS[module]
        command = ["node", str(ROOT / test.check), "--package", str(wasm32_package(module))]
        if test.wasm64:
            command += ["--package", str(ROOT / "www" / "public" / "pkg_wasm64" / module)]
        for other in test.uses:
            command += [f"--{other}-package", str(wasm32_package(other))]
        command += ["--update-golden"] if args.update_golden else []
        command += [a for a in args.check_args if a != "--"]
        if subprocess.run(command).returncode != 0:
            failed.append(module)
    print(f"WASM checks failed: {', '.join(failed)}" if failed else "All WASM checks passed.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
