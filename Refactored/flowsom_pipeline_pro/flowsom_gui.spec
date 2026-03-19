# flowsom_gui.spec — PyInstaller spec optimisé pour FlowSOM Analyzer Pro (GUI)
#
# Utilisation :
#   cd "C:\Users\Florian Travail\Documents\FlowSom\Perplexity\flowsom_pipeline_pro"
#   python -m PyInstaller flowsom_gui.spec -y
#
# Mode onedir : résultat dans dist/FlowSOMAnalyzer/FlowSOMAnalyzer.exe
# → lancement instantané (pas d'extraction au démarrage)
# → pour distribuer, zipper le dossier dist/FlowSOMAnalyzer/ entier
#
# OPTIMISATIONS v4 :
#   - collect_all() pour les bibliothèques scientifiques complexes
#     (anndata, zarr, scanpy, scipy, sklearn, numpy, pandas, igraph, numba,
#      umap-learn, networkx, h5py, matplotlib, plotly, flowsom, readfcs, mudata)
#   - Exclut torch (~4 GB) et son écosystème
#   - Inclut rich (requis par zarr.core._tree)
#   - Inclut leidenalg (requis par scanpy.tools._leiden)
#   - UPX désactivé (non installé)

import sys
from pathlib import Path
from PyInstaller.utils.hooks import collect_data_files, collect_submodules, copy_metadata, collect_all

_HERE = Path(SPECPATH)  # répertoire du .spec (= racine du package)

# ── Collect-all pour les bibliothèques complexes ───────────────────────────
# collect_all() retourne (datas, binaries, hiddenimports) — beaucoup plus
# fiable que collect_submodules() seul pour les gros packages scientifiques.
datas = []
binaries = []
hidden_imports = []

for _pkg in (
    "anndata", "zarr", "scanpy", "scipy", "sklearn", "numpy", "pandas",
    "networkx", "igraph", "numba", "umap", "plotly", "matplotlib",
    "flowsom", "readfcs", "mudata", "flowio", "h5py", "loguru",
    "seaborn", "tqdm", "PIL", "yaml", "natsort", "packaging",
    "fsspec", "rich",
    # Monitoring de performance
    "psutil", "GPUtil",
    # PDF + export statique Plotly
    "reportlab", "kaleido",
):
    try:
        _d, _b, _h = collect_all(_pkg)
        datas += _d
        binaries += _b
        hidden_imports += _h
    except Exception:
        pass

# ── Données à embarquer ────────────────────────────────────────────────────
datas += [
    # YAML de configuration par défaut → déposé dans config/ à côté de l'exe
    (str(_HERE / "config" / "default_config.yaml"), "config"),
    # Configuration MRD résiduelle
    (str(_HERE / "config" / "mrd_config.yaml"), "config"),
]

# Metadata package pour importlib.metadata.version(...)
# Beaucoup de packages (scanpy, flowsom, anndata…) appellent
# importlib.metadata.version("X") au démarrage → il faut embarquer le dist-info.
for _pkg in (
    "flowsom", "anndata", "scanpy", "scikit-learn", "scipy", "numpy",
    "pandas", "matplotlib", "seaborn", "plotly", "umap-learn", "numba",
    "igraph", "networkx", "h5py", "loguru", "mudata", "readfcs",
    "flowio", "fcswrite", "PyYAML", "tqdm", "Pillow", "kaleido",
    "natsort", "packaging", "zarr", "fsspec", "rich", "leidenalg",
    # Monitoring + PDF
    "psutil", "GPUtil", "reportlab",
):
    try:
        datas += copy_metadata(_pkg)
    except Exception:
        pass

# PyQtWebEngine : ressources Chromium (nécessaires pour QWebEngineView)
try:
    datas += collect_data_files("PyQt5.QtWebEngineWidgets")
except Exception:
    pass
try:
    datas += collect_data_files("PyQtWebEngine")
except Exception:
    pass

# ── Hidden imports ─────────────────────────────────────────────────────────
hidden_imports += [
    # ── PyQt5 + WebEngine ──────────────────────────────────────────────────
    "PyQt5",
    "PyQt5.QtCore",
    "PyQt5.QtGui",
    "PyQt5.QtWidgets",
    "PyQt5.QtSvg",
    "PyQt5.QtPrintSupport",
    "PyQt5.QtWebEngineWidgets",
    "PyQt5.QtWebEngineCore",
    "PyQt5.QtWebChannel",
    "PyQt5.QtNetwork",
    # ── Matplotlib backends ────────────────────────────────────────────────
    "matplotlib.backends.backend_qt5agg",
    "matplotlib.backends.backend_agg",
    "PIL.ImageQt",
    "PIL.ImageDraw",
    "mpl_toolkits",
    # ── NumPy internals ────────────────────────────────────────────────────
    "numpy.core._multiarray_umath",
    "numpy.core.multiarray",
    # ── SciPy internals ───────────────────────────────────────────────────
    "scipy.spatial.distance",
    "scipy.cluster.hierarchy",
    "scipy.sparse.csgraph",
    "scipy._lib.messagestream",
    # ── Scikit-learn internals ────────────────────────────────────────────
    "sklearn.utils._cython_blas",
    "sklearn.neighbors._partition_nodes",
    "sklearn.utils._typedefs",
    "sklearn.utils._heap",
    "sklearn.utils._sorting",
    "sklearn.utils._vector_sentinel",
    "sklearn.tree._utils",
    # ── Pandas internals ──────────────────────────────────────────────────
    "pandas._libs.tslibs.timedeltas",
    "pandas._libs.tslibs.np_datetime",
    "pandas._libs.tslibs.nattype",
    "pandas._libs.tslibs.base",
    "pandas._libs.skiplist",
    # ── AnnData ───────────────────────────────────────────────────────────
    "anndata._core.merge",
    # ── FlowSOM / FlowIO / readfcs ────────────────────────────────────────
    "readfcs._core",
    "readfcs.datasets",
    # ── Scanpy (requis par flowsom.pl.plot_functions au top-level) ────────
    "scanpy.preprocessing",
    "scanpy.tools",
    # ── Leiden (requis par scanpy.tools._leiden) ──────────────────────────
    "leidenalg",
    # ── Rich (requis par zarr.core._tree) ────────────────────────────────
    "rich.console",
    "rich.tree",
    # ── UMAP + numba ──────────────────────────────────────────────────────
    "llvmlite",
    "llvmlite.binding",
    # ── Monitoring de performance ─────────────────────────────────────────
    "psutil",
    "psutil._pswindows",
    "psutil._pslinux",
    "psutil._psosx",
    "GPUtil",
    # ── PDF (ReportLab) — imports internes non détectés par PyInstaller ───
    "reportlab",
    "reportlab.pdfgen",
    "reportlab.pdfgen.canvas",
    "reportlab.lib",
    "reportlab.lib.pagesizes",
    "reportlab.lib.styles",
    "reportlab.lib.units",
    "reportlab.lib.colors",
    "reportlab.lib.enums",
    "reportlab.lib.utils",
    "reportlab.platypus",
    "reportlab.platypus.flowables",
    "reportlab.platypus.frames",
    "reportlab.platypus.paragraph",
    "reportlab.platypus.tables",
    "reportlab.platypus.doctemplate",
    "reportlab.graphics",
    "reportlab.rl_config",
    # ── Optionnels ────────────────────────────────────────────────────────
    "fcswrite",
    "fcsparser",
    "pytometry",
    "certifi",
    "pynndescent",
]

# Collecte automatique du projet lui-même
hidden_imports += collect_submodules("flowsom_pipeline_pro")

# ── Analyse ───────────────────────────────────────────────────────────────
a = Analysis(
    [str(_HERE / "launch_gui.py")],
    pathex=[str(_HERE), str(_HERE.parent)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[str(_HERE / "runtime_hook_debug.py")],
    excludes=[
        # ── UI alternatives inutiles ──────────────────────────────────────
        "tkinter",
        "wx",
        "gi",
        # ── Jupyter / notebooks (inutiles) ────────────────────────────────
        "IPython",
        "jupyter",
        "notebook",
        "nbformat",
        "nbconvert",
        "ipywidgets",
        "tornado",
        "zmq",
        # ── TORCH et tout son écosystème (~4 GB !) ────────────────────────
        # torch n'est PAS utilisé par flowsom_pipeline_pro.
        # Scanpy et umap l'importent de façon OPTIONNELLE (try/except).
        "torch",
        "torchvision",
        "torchaudio",
        "torch.distributed",
        "torch.cuda",
        "torch.nn",
        "torch.optim",
        "torch.utils",
        "torch.onnx",
        "functorch",
        # ── ONNX / ONNXRuntime (inutiles) ─────────────────────────────────
        "onnxruntime",
        "onnx",
        "onnxscript",
        # ── Selenium / web scraping (inutile) ─────────────────────────────
        "selenium",
        # ── Tests (inutiles en prod) ───────────────────────────────────────
        "hypothesis",
        "pytest",
        "test",
        "tests",
        # ── Visualisation lourde non utilisée ─────────────────────────────
        "bokeh",
        # ── Analytics lourds non utilisés directement ─────────────────────
        "statsmodels",
        "pyarrow",
        # ── Autres non utilisés ───────────────────────────────────────────
        "dask",
        "xarray",
        "sympy",
        "cvxpy",
        "numexpr",
        "tables",
    ],
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data)

# ── Mode onedir : EXE léger + dossier _internal — démarrage instantané ────
exe = EXE(
    pyz,
    a.scripts,
    [],                   # Binaires/datas gérés par COLLECT (pas onefile)
    exclude_binaries=True,
    name="FlowSOMAnalyzer",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,            # UPX désactivé (non installé)
    console=False,        # Pas de fenêtre console (GUI pur)
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    name="FlowSOMAnalyzer",   # → dist/FlowSOMAnalyzer/FlowSOMAnalyzer.exe
)
