"""
launch_gui.py — Point d'entrée pour la compilation PyInstaller (GUI uniquement).

Usage:
    python launch_gui.py       (développement)
    pyinstaller flowsom_gui.spec  (production → .exe)
"""

import sys
import os
import logging as _logging
from pathlib import Path


def _enable_windows_crisp_rendering() -> None:
    """Request per-monitor DPI awareness to avoid blurry bitmap-scaled UI on Windows."""
    if os.name != "nt":
        return

    # Qt DPI env hints must be set before QApplication is created.
    os.environ.setdefault("QT_AUTO_SCREEN_SCALE_FACTOR", "1")
    os.environ.setdefault("QT_ENABLE_HIGHDPI_SCALING", "1")
    os.environ.setdefault("QT_SCALE_FACTOR_ROUNDING_POLICY", "RoundPreferFloor")

    try:
        import ctypes

        # PER_MONITOR_AWARE_V2 gives the sharpest text on mixed-DPI displays.
        ctypes.windll.user32.SetProcessDpiAwarenessContext(ctypes.c_void_p(-4))
        return
    except Exception:
        pass

    try:
        import ctypes

        ctypes.windll.shcore.SetProcessDpiAwareness(2)
        return
    except Exception:
        pass

    try:
        import ctypes

        ctypes.windll.user32.SetProcessDPIAware()
    except Exception:
        pass


# ── Mode windowed (console=False) : sys.stderr/stdout sont None → crash logging ─
# On redirige vers un flux nul AVANT tout import pour éviter les
# AttributeError: 'NoneType' object has no attribute 'write'.
if getattr(sys, "frozen", False):
    if sys.stderr is None:
        sys.stderr = open(os.devnull, "w", encoding="utf-8")
    if sys.stdout is None:
        sys.stdout = open(os.devnull, "w", encoding="utf-8")

# ── Résolution du chemin vers les dépendances dans le .exe ─────────────────────
# En mode onedir, sys._MEIPASS = sous-dossier _internal à côté de l'exe.
if getattr(sys, "frozen", False):
    _BASE_DIR = Path(sys._MEIPASS)
else:
    _BASE_DIR = Path(__file__).resolve().parent

# Rendre le package flowsom_pipeline_pro importable (mode frozen ou dev)
if str(_BASE_DIR) not in sys.path:
    sys.path.insert(0, str(_BASE_DIR))

# Idem pour le dossier parent (mode dev sans pip install)
_PARENT = _BASE_DIR.parent
if str(_PARENT) not in sys.path:
    sys.path.insert(0, str(_PARENT))

# ── Suppression des logs verbeux de kaleido (Chromium/CDP) ─────────────────────
# kaleido utilise le logger root Python — on le filtre à WARNING pour ne garder
# que les vraies erreurs, pas les traces internes Chromium/CDP.
for _noisy in (
    "kaleido",
    "kaleido.scopes",
    "kaleido.scopes.base",
    "kaleido.scopes.chromium",
    "chromote",
    "pyppeteer",
    "matplotlib",
    "matplotlib.font_manager",
):
    _logging.getLogger(_noisy).setLevel(_logging.WARNING)

# Important: must run before importing Qt modules / creating QApplication.
_enable_windows_crisp_rendering()

from flowsom_pipeline_pro.gui.main_window import main

if __name__ == "__main__":
    main()
