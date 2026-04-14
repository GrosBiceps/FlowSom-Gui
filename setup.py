"""
setup.py — Packaging de flowsom_pipeline_pro.

Installation:
    pip install -e .                      # mode développement (CPU, sans GUI)
    pip install -e ".[gui]"               # avec interface PyQt5
    pip install -e ".[gpu]"               # avec PyTorch CUDA + scanpy
    pip install -e ".[full]"              # toutes les dépendances CPU
    pip install -e ".[full,gpu]"          # installation complète (CPU + GPU)

Mode dégradé :
    L'application fonctionne sans torch (GPU désactivé) et sans scanpy.
    Les fonctionnalités GPU sont automatiquement ignorées si torch est absent.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Lecture du README pour la description longue
_HERE = Path(__file__).parent
_README = (
    (_HERE / "README.md").read_text(encoding="utf-8")
    if (_HERE / "README.md").exists()
    else ""
)

setup(
    name="flowsom_pipeline_pro",
    version="1.0.0",
    description=(
        "Pipeline d'analyse de cytométrie en flux MRD (Maladie Résiduelle Détectable) "
        "basé sur FlowSOM — architecture modulaire production-ready"
    ),
    long_description=_README,
    long_description_content_type="text/markdown",
    author="Florian",
    python_requires=">=3.10",
    # Découverte automatique de tous les sous-paquets
    packages=find_packages(exclude=["tests*", "*.tests*"]),
    # Inclure les fichiers de données (YAML par défaut)
    package_data={
        "flowsom_pipeline_pro": [
            "config/default_config.yaml",
            "config/mrd_config.yaml",
            "config/panels/*.yaml",
        ],
    },
    # ── Dépendances obligatoires ───────────────────────────────────────────────
    install_requires=[
        "numpy>=1.24",
        "pandas>=2.0",
        "scikit-learn>=1.3",
        "scipy>=1.11",
        "matplotlib>=3.7",
        "seaborn>=0.12",
        "anndata>=0.10",
        "pyyaml>=6.0",
        "tqdm>=4.65",
    ],
    # ── Dépendances optionnelles ───────────────────────────────────────────────
    extras_require={
        # FlowSOM principal (saeyslab)
        "flowsom": ["flowsom==0.2.2"],
        # Lecture/écriture FCS
        "fcs": ["flowio==1.4.0"],
        # Export FCS avec colonnes de clustering
        "fcs_export": ["fcswrite==0.6.4"],
        # Transformation logicle via FlowKit
        "logicle": ["flowkit"],
        # Transformations cytométriques (pytometry = scanpy-based)
        "pytometry": ["pytometry==0.1.6"],
        # Interface graphique PyQt5
        "gui": [
            "PyQt5==5.15.11",
            "PyQtWebEngine==5.15.7",
            "qtawesome>=1.3.0",
        ],
        # Accélération GPU — PyTorch CUDA (choisir la variante CUDA via --index-url)
        # Ne pas installer automatiquement : nécessite une commande pip spécifique.
        # Voir requirements-gpu.txt pour les instructions détaillées.
        "gpu": [
            "torch>=2.0",
            "cupy",
            "scanpy==1.11.5",
        ],
        # Notebooks Jupyter (développement / exploration)
        "notebooks": [
            "ipykernel>=6.0",
            "ipywidgets>=8.0",
            "tqdm>=4.65",
        ],
        # Génération de rapports PDF
        "reports": ["reportlab>=4.0"],
        # Toutes les dépendances CPU (sans GPU, sans notebooks)
        "full": [
            "flowsom==0.2.2",
            "flowio==1.4.0",
            "fcswrite==0.6.4",
            "pytometry==0.1.6",
            "umap-learn>=0.5",
            "anndata>=0.10",
            "mudata>=0.3.0",
            "igraph>=1.0.0",
            "networkx>=3.0",
            "reportlab>=4.0",
            "kaleido>=0.1.0",
            "plotly>=5.0",
            "PyQt5==5.15.11",
            "PyQtWebEngine==5.15.7",
            "qtawesome>=1.3.0",
            "psutil>=5.9",
            "GPUtil>=1.4.0",
            "loguru>=0.7",
        ],
    },
    # ── Points d'entrée CLI ───────────────────────────────────────────────────
    entry_points={
        "console_scripts": [
            "flowsom-analyze=flowsom_pipeline_pro.cli.main:main",
        ],
    },
    # ── Classifieurs PyPI ─────────────────────────────────────────────────────
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    keywords=[
        "flowsom",
        "cytometry",
        "flow cytometry",
        "MRD",
        "AML",
        "bioinformatics",
        "clustering",
    ],
)
