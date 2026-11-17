"""Compile (if needed) and load the NMODL mechanisms in neuron/.

Run from the project venv (Python 3.13, pip `neuron==9.0.2`) with PYTHONPATH unset: the app-bundle NEURON in
PYTHONPATH shadows pip's and cannot load under Homebrew's Python (T0, RESULTS.md section 1).
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MOD_DIR = ROOT / "neuron"
_ARCHS = ("arm64", "x86_64", "aarch64")
_loaded = False


def _find_lib():
    for arch in _ARCHS:
        for sub in ("", ".libs"):
            for ext in ("dylib", "so"):
                lib = MOD_DIR / arch / sub / f"libnrnmech.{ext}"
                if lib.exists():
                    return lib
    return None


def load_mechanisms():
    """Run nrnivmodl on neuron/ when the build is missing or older than a .mod file, then load it once."""
    global _loaded
    from neuron import h

    if _loaded:
        return h
    lib = _find_lib()
    newest_mod = max(p.stat().st_mtime for p in MOD_DIR.glob("*.mod"))
    if lib is None or lib.stat().st_mtime < newest_mod:
        nrnivmodl = Path(sys.executable).parent / "nrnivmodl"
        subprocess.run([str(nrnivmodl) if nrnivmodl.exists() else "nrnivmodl"],
                       cwd=MOD_DIR, check=True, capture_output=True)
        lib = _find_lib()
    if lib is None or not h.nrn_load_dll(str(lib)):
        raise RuntimeError(f"could not load NMODL library from {MOD_DIR}")
    h.load_file("stdrun.hoc")
    _loaded = True
    return h
