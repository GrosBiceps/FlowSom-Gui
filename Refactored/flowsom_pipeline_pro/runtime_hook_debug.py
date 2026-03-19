"""
Runtime hook — Redirection stderr/stdout + diagnostic d'imports pour le mode frozen.
S'exécute AVANT le script principal (launch_gui.py).
"""
import sys
import os
import traceback
from pathlib import Path

# ── Redirection stderr/stdout IMMÉDIATE (avant tout import de logging) ────────
# En mode windowed (console=False), sys.stderr/stdout sont None → crash des
# StreamHandlers de logging. On redirige vers devnull au plus tôt.
if getattr(sys, "frozen", False):
    if sys.stderr is None:
        sys.stderr = open(os.devnull, "w", encoding="utf-8")
    if sys.stdout is None:
        sys.stdout = open(os.devnull, "w", encoding="utf-8")

def _debug_imports():
    if not getattr(sys, "frozen", False):
        return

    # Log à côté de l'exe
    exe_dir = Path(sys.executable).parent
    log_path = exe_dir / "import_debug.log"

    lines = []
    lines.append(f"=== Import Debug — {__import__('datetime').datetime.now()} ===")
    lines.append(f"sys.executable: {sys.executable}")
    lines.append(f"sys._MEIPASS:   {getattr(sys, '_MEIPASS', 'N/A')}")
    lines.append(f"sys.path:       {sys.path[:5]}")
    lines.append("")

    # Test chaque dépendance critique de flowsom dans l'ordre
    test_chain = [
        "numpy",
        "pandas",
        "scipy",
        "anndata",
        "flowio",
        "readfcs",
        "readfcs._core",
        "mudata",
        "loguru",
        "igraph",
        "networkx",
        "scanpy",
        "scanpy.preprocessing",
        "scanpy.tools",
        "numba",
        "flowsom",
        "flowsom.io",
        "flowsom.models",
        "flowsom.pl",
        "flowsom.pp",
        "flowsom.tl",
        "flowsom.main",
    ]

    for mod_name in test_chain:
        try:
            __import__(mod_name)
            lines.append(f"OK     {mod_name}")
        except Exception as exc:
            lines.append(f"FAIL   {mod_name}")
            lines.append(f"       {type(exc).__name__}: {exc}")
            lines.append(f"       {''.join(traceback.format_tb(exc.__traceback__)[-2:]).strip()}")
            lines.append("")

    lines.append("\n=== Fin du diagnostic ===")

    try:
        log_path.write_text("\n".join(lines), encoding="utf-8")
    except Exception:
        pass

_debug_imports()
