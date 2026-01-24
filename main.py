# -*- coding: utf-8 -*-

import sys
import os
import tempfile
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
import warnings
warnings.filterwarnings('ignore')

# Imports PyQt5 pour l'interface graphique
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QFileDialog, QGroupBox, QSpinBox,
    QDoubleSpinBox, QCheckBox, QComboBox, QProgressBar, QStatusBar,
    QMessageBox, QTabWidget, QScrollArea, QFrame, QSplitter,
    QTableWidget, QTableWidgetItem, QHeaderView, QGridLayout,
    QLineEdit, QTextEdit, QSizePolicy, QGraphicsDropShadowEffect,
    QListWidget, QListWidgetItem, QAbstractItemView
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QSize, QPropertyAnimation, QEasingCurve
from PyQt5.QtGui import QFont, QColor, QPalette, QIcon, QLinearGradient, QPainter, QBrush

# Imports pour le graphique Matplotlib intégré
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Imports scientifiques
import numpy as np
import pandas as pd
import anndata as ad

# Import FlowSOM (saeyslab)
try:
    import flowsom as fs
    FLOWSOM_AVAILABLE = True
except ImportError:
    FLOWSOM_AVAILABLE = False
    print("ATTENTION: La librairie flowsom n'est pas installee.")
    print("   Installez-la avec: pip install flowsom")

# Import pour l'export FCS
try:
    import fcswrite
    FCSWRITE_AVAILABLE = True
except ImportError:
    FCSWRITE_AVAILABLE = False
    print("ATTENTION: La librairie fcswrite n'est pas installee.")
    print("   Installez-la avec: pip install fcswrite")

# Import pour les métriques de clustering
from sklearn.metrics import silhouette_score
from sklearn.cluster import AgglomerativeClustering

# Import scanpy pour UMAP/t-SNE
try:
    import scanpy as sc
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False
    print("ATTENTION: La librairie scanpy n'est pas installee.")
    print("   Installez-la avec: pip install scanpy")

# Import UMAP directement (plus stable que via scanpy)
try:
    import umap
    UMAP_AVAILABLE = True
except ImportError:
    UMAP_AVAILABLE = False

# Import sklearn pour t-SNE (plus stable)
try:
    from sklearn.manifold import TSNE
    TSNE_AVAILABLE = True
except ImportError:
    TSNE_AVAILABLE = False

# Import FlowKit pour les transformations Logicle
try:
    import flowkit as fk
    FLOWKIT_AVAILABLE = True
except ImportError:
    FLOWKIT_AVAILABLE = False
    print("INFO: FlowKit non installe (transformations Logicle).")
    print("   pip install flowkit")

# Import ReportLab pour la generation de PDF
try:
    from reportlab.lib.pagesizes import A4, letter
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import cm, inch
    from reportlab.lib.colors import HexColor
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
    from reportlab.lib import colors
    REPORTLAB_AVAILABLE = True
except ImportError:
    REPORTLAB_AVAILABLE = False
    print("INFO: ReportLab non installe (export PDF).")
    print("   pip install reportlab")

# Import scipy pour les statistiques
from scipy import stats

# Import datetime pour les horodatages
from datetime import datetime
import json

# Import YAML pour la configuration
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False
    print("INFO: PyYAML non installe (configuration YAML).")
    print("   pip install pyyaml")


# =============================================================================
# CHARGEMENT DE LA CONFIGURATION
# =============================================================================

def load_config() -> Dict[str, Any]:
    """
    Charge la configuration depuis config.yaml.
    Si le fichier n'existe pas ou YAML non disponible, retourne la config par defaut.
    """
    default_config = {
        'modules_principaux': {
            'flowsom_analysis': True,
            'visualization_tab': True,
            'statistics_tab': True,
            'export_functions': True
        },
        'onglets': {
            'clusters_interactifs': True,
            'fcs_visualization': True,
            'suivi_patient': True,
            'logs': True
        },
        'modules_optionnels': {
            'data_transformations': True,
            'patient_tracking': True,
            'umap_tsne': True,
            'interactive_clusters': True,
            'flowsom_mst_stars': True,
            'pdf_reports': True
        },
        'modules_avances': {
            'nbm_reference': True,
            'nbm_comparison': True,
            'clone_persistence': True
        },
        'ui_settings': {
            'theme': 'dark',
            'font_size': 10,
            'show_tooltips': True,
            'default_embedding_cells': 5000,
            'show_references': True
        }
    }
    
    if not YAML_AVAILABLE:
        return default_config
    
    config_path = Path(__file__).parent / 'config.yaml'
    if not config_path.exists():
        return default_config
    
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            user_config = yaml.safe_load(f)
        
        # Merge avec la config par defaut
        def deep_merge(base: Dict, overlay: Dict) -> Dict:
            result = base.copy()
            for key, value in overlay.items():
                if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                    result[key] = deep_merge(result[key], value)
                else:
                    result[key] = value
            return result
        
        return deep_merge(default_config, user_config)
    except Exception as e:
        print(f"Erreur chargement config.yaml: {e}")
        return default_config


# Configuration globale
CONFIG = load_config()


# =============================================================================
# STYLES CSS MODERNES POUR L'INTERFACE
# =============================================================================

STYLESHEET = """
/* ============================================
   THÈME PRINCIPAL - CATPPUCCIN MOCHA AMÉLIORÉ
   ============================================ */

QMainWindow {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 #11111b, stop:0.5 #1e1e2e, stop:1 #181825);
}

QWidget {
    background-color: transparent;
    color: #cdd6f4;
    font-family: 'Segoe UI', 'SF Pro Display', Arial, sans-serif;
    font-size: 10pt;
}

QWidget#centralWidget {
    background: transparent;
}

QWidget#leftPanel {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 rgba(49, 50, 68, 0.95), stop:1 rgba(30, 30, 46, 0.95));
    border-right: 1px solid rgba(137, 180, 250, 0.3);
    border-radius: 0px;
}

QWidget#rightPanel {
    background: rgba(30, 30, 46, 0.7);
}

/* ============================================
   GROUPBOX AVEC EFFET GLASSMORPHISM
   ============================================ */

QGroupBox {
    font-weight: 600;
    font-size: 11pt;
    border: 1px solid rgba(137, 180, 250, 0.2);
    border-radius: 12px;
    margin-top: 20px;
    padding: 20px 15px 15px 15px;
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 rgba(69, 71, 90, 0.6), stop:1 rgba(49, 50, 68, 0.8));
}

QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 15px;
    top: 5px;
    padding: 5px 15px;
    color: #89b4fa;
    background: linear-gradient(90deg, rgba(137, 180, 250, 0.2), transparent);
    border-radius: 8px;
    font-size: 11pt;
}

/* ============================================
   BOUTONS AVEC ANIMATIONS ET GRADIENTS
   ============================================ */

QPushButton {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #585b70, stop:1 #45475a);
    border: 1px solid rgba(255, 255, 255, 0.1);
    border-radius: 10px;
    padding: 12px 24px;
    color: #cdd6f4;
    font-weight: 600;
    font-size: 10pt;
    min-height: 20px;
}

QPushButton:hover {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #6c7086, stop:1 #585b70);
    border: 1px solid rgba(137, 180, 250, 0.4);
}

QPushButton:pressed {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #45475a, stop:1 #313244);
}

QPushButton:disabled {
    background: rgba(49, 50, 68, 0.5);
    color: #6c7086;
    border: 1px solid rgba(255, 255, 255, 0.05);
}

QPushButton#primaryBtn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 #89b4fa, stop:1 #74c7ec);
    color: #11111b;
    border: none;
    font-weight: 700;
}

QPushButton#primaryBtn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 #b4befe, stop:1 #89dceb);
}

QPushButton#primaryBtn:disabled {
    background: rgba(137, 180, 250, 0.3);
    color: rgba(17, 17, 27, 0.5);
}

QPushButton#successBtn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 #a6e3a1, stop:1 #94e2d5);
    color: #11111b;
    border: none;
    font-weight: 600;
}

QPushButton#successBtn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 #b4f5b0, stop:1 #a6f3e8);
}

QPushButton#dangerBtn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 #f38ba8, stop:1 #fab387);
    color: #11111b;
    border: none;
    font-weight: 600;
}

QPushButton#dangerBtn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 #f5a0ba, stop:1 #fcc8a0);
}

QPushButton#exportBtn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 #cba6f7, stop:1 #f5c2e7);
    color: #11111b;
    border: none;
    font-weight: 600;
}

QPushButton#exportBtn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 #dbb8ff, stop:1 #f8d5ed);
}

/* ============================================
   LABELS
   ============================================ */

QLabel {
    color: #cdd6f4;
    background: transparent;
}

QLabel#titleLabel {
    font-size: 22pt;
    font-weight: 700;
    color: transparent;
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 #89b4fa, stop:0.5 #cba6f7, stop:1 #f5c2e7);
    background-clip: text;
    -webkit-background-clip: text;
}

QLabel#subtitleLabel {
    font-size: 10pt;
    color: #a6adc8;
    font-style: italic;
}

QLabel#sectionLabel {
    font-size: 11pt;
    font-weight: 600;
    color: #89b4fa;
}

QLabel#fileLabel {
    font-size: 9pt;
    color: #a6adc8;
    padding: 5px 10px;
    background: rgba(69, 71, 90, 0.4);
    border-radius: 6px;
}

QLabel#statusSuccess {
    color: #a6e3a1;
    font-weight: 600;
}

QLabel#statusError {
    color: #f38ba8;
    font-weight: 600;
}

/* ============================================
   INPUTS (SPINBOX, COMBOBOX)
   ============================================ */

QSpinBox, QDoubleSpinBox, QComboBox, QLineEdit {
    background: rgba(69, 71, 90, 0.6);
    border: 1px solid rgba(137, 180, 250, 0.2);
    border-radius: 8px;
    padding: 10px 12px;
    color: #cdd6f4;
    min-height: 20px;
    font-size: 10pt;
    selection-background-color: #89b4fa;
}

QSpinBox:hover, QDoubleSpinBox:hover, QComboBox:hover, QLineEdit:hover {
    border: 1px solid rgba(137, 180, 250, 0.4);
    background: rgba(69, 71, 90, 0.8);
}

QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus, QLineEdit:focus {
    border: 2px solid #89b4fa;
    background: rgba(69, 71, 90, 0.9);
}

QSpinBox::up-button, QSpinBox::down-button,
QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {
    background: rgba(137, 180, 250, 0.2);
    border: none;
    border-radius: 4px;
    width: 20px;
    margin: 2px;
}

QSpinBox::up-button:hover, QSpinBox::down-button:hover,
QDoubleSpinBox::up-button:hover, QDoubleSpinBox::down-button:hover {
    background: rgba(137, 180, 250, 0.4);
}

QComboBox::drop-down {
    border: none;
    width: 30px;
    background: rgba(137, 180, 250, 0.2);
    border-top-right-radius: 8px;
    border-bottom-right-radius: 8px;
}

QComboBox::down-arrow {
    width: 12px;
    height: 12px;
}

QComboBox QAbstractItemView {
    background: #313244;
    border: 1px solid rgba(137, 180, 250, 0.3);
    border-radius: 8px;
    selection-background-color: #89b4fa;
    selection-color: #11111b;
    padding: 5px;
}

/* ============================================
   CHECKBOX MODERNE
   ============================================ */

QCheckBox {
    spacing: 10px;
    color: #cdd6f4;
    font-size: 10pt;
}

QCheckBox::indicator {
    width: 22px;
    height: 22px;
    border-radius: 6px;
    border: 2px solid rgba(137, 180, 250, 0.3);
    background: rgba(69, 71, 90, 0.4);
}

QCheckBox::indicator:hover {
    border-color: rgba(137, 180, 250, 0.6);
    background: rgba(69, 71, 90, 0.6);
}

QCheckBox::indicator:checked {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 #89b4fa, stop:1 #74c7ec);
    border-color: #89b4fa;
}

/* ============================================
   PROGRESSBAR ANIMÉE
   ============================================ */

QProgressBar {
    border: none;
    border-radius: 10px;
    background: rgba(69, 71, 90, 0.4);
    text-align: center;
    color: #cdd6f4;
    font-weight: 600;
    font-size: 9pt;
    min-height: 20px;
}

QProgressBar::chunk {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 #89b4fa, stop:0.5 #cba6f7, stop:1 #f5c2e7);
    border-radius: 10px;
}

/* ============================================
   STATUS BAR
   ============================================ */

QStatusBar {
    background: rgba(49, 50, 68, 0.9);
    color: #a6adc8;
    border-top: 1px solid rgba(137, 180, 250, 0.2);
    font-size: 9pt;
    padding: 5px;
}

/* ============================================
   TABWIDGET MODERNE
   ============================================ */

QTabWidget::pane {
    border: 1px solid rgba(137, 180, 250, 0.2);
    border-radius: 12px;
    background: rgba(49, 50, 68, 0.6);
    margin-top: 0px;
    padding-top: 5px;
}

QTabWidget::tab-bar {
    left: 10px;
}

QTabBar {
    qproperty-drawBase: 0;
}

QTabBar::tab {
    background: rgba(69, 71, 90, 0.6);
    border: 1px solid rgba(137, 180, 250, 0.15);
    border-bottom: none;
    padding: 10px 18px;
    margin-right: 4px;
    margin-top: 3px;
    border-top-left-radius: 8px;
    border-top-right-radius: 8px;
    color: #a6adc8;
    font-weight: 500;
    font-size: 9pt;
    min-width: 80px;
}

QTabBar::tab:selected {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 #89b4fa, stop:1 #cba6f7);
    color: #11111b;
    font-weight: 600;
    margin-top: 0px;
    padding-bottom: 12px;
}

QTabBar::tab:hover:!selected {
    background: rgba(137, 180, 250, 0.3);
    color: #cdd6f4;
}

QTabBar::tab:first {
    margin-left: 0px;
}

/* ============================================
   SCROLLAREA & SCROLLBAR
   ============================================ */

QScrollArea {
    border: none;
    background: transparent;
}

QScrollBar:vertical {
    background: rgba(69, 71, 90, 0.3);
    width: 10px;
    border-radius: 5px;
    margin: 0;
}

QScrollBar::handle:vertical {
    background: rgba(137, 180, 250, 0.4);
    border-radius: 5px;
    min-height: 30px;
}

QScrollBar::handle:vertical:hover {
    background: rgba(137, 180, 250, 0.6);
}

QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0;
}

QScrollBar:horizontal {
    background: rgba(69, 71, 90, 0.3);
    height: 10px;
    border-radius: 5px;
}

QScrollBar::handle:horizontal {
    background: rgba(137, 180, 250, 0.4);
    border-radius: 5px;
    min-width: 30px;
}

QScrollBar::handle:horizontal:hover {
    background: rgba(137, 180, 250, 0.6);
}

/* ============================================
   TABLEWIDGET
   ============================================ */

QTableWidget {
    background: rgba(49, 50, 68, 0.6);
    border: 1px solid rgba(137, 180, 250, 0.2);
    border-radius: 10px;
    gridline-color: rgba(137, 180, 250, 0.1);
    font-size: 10pt;
}

QTableWidget::item {
    padding: 10px;
    color: #cdd6f4;
    border-bottom: 1px solid rgba(137, 180, 250, 0.1);
}

QTableWidget::item:selected {
    background: rgba(137, 180, 250, 0.3);
    color: #cdd6f4;
}

QTableWidget::item:hover {
    background: rgba(137, 180, 250, 0.15);
}

QHeaderView::section {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 rgba(137, 180, 250, 0.3), stop:1 rgba(137, 180, 250, 0.1));
    color: #cdd6f4;
    padding: 12px;
    border: none;
    font-weight: 600;
    font-size: 10pt;
}

/* ============================================
   TEXTEDIT (LOGS)
   ============================================ */

QTextEdit {
    background: rgba(17, 17, 27, 0.8);
    border: 1px solid rgba(137, 180, 250, 0.2);
    border-radius: 10px;
    color: #a6e3a1;
    padding: 12px;
    font-family: 'Consolas', 'Fira Code', monospace;
    font-size: 9pt;
    selection-background-color: #89b4fa;
}

/* ============================================
   SPLITTER
   ============================================ */

QSplitter::handle {
    background: rgba(137, 180, 250, 0.2);
    width: 3px;
}

QSplitter::handle:hover {
    background: rgba(137, 180, 250, 0.5);
}

/* ============================================
   TOOLTIPS
   ============================================ */

QToolTip {
    background: #313244;
    color: #cdd6f4;
    border: 1px solid rgba(137, 180, 250, 0.3);
    border-radius: 6px;
    padding: 8px 12px;
    font-size: 9pt;
}
"""


# =============================================================================
# CLASSES UTILITAIRES POUR LE PROTOCOLE LAM MRD
# =============================================================================

class PatientData:
    """
    Gestion des donnees patient avec suivi longitudinal (Time-Series).
    """
    
    def __init__(self, patient_id: str):
        self.patient_id = patient_id
        self.timepoints: Dict[str, Dict[str, Any]] = {}  # {timepoint_id: data}
        self.metadata: Dict[str, Any] = {}
        self.created_at = datetime.now().isoformat()
        
    def add_timepoint(self, timepoint_id: str, files: List[str], 
                      date: Optional[str] = None, notes: str = ""):
        """Ajoute un point de temps avec ses fichiers FCS."""
        self.timepoints[timepoint_id] = {
            'files': files,
            'date': date or datetime.now().strftime("%Y-%m-%d"),
            'notes': notes,
            'results': None,
            'lsc_score': None
        }
        
    def set_timepoint_results(self, timepoint_id: str, results: Dict[str, Any]):
        """Enregistre les resultats d'analyse pour un point de temps."""
        if timepoint_id in self.timepoints:
            self.timepoints[timepoint_id]['results'] = results
            
    def get_timeline(self) -> List[Tuple[str, str]]:
        """Retourne la timeline ordonnee par date."""
        items = [(tid, tp['date']) for tid, tp in self.timepoints.items()]
        return sorted(items, key=lambda x: x[1])
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialise les donnees patient."""
        return {
            'patient_id': self.patient_id,
            'timepoints': self.timepoints,
            'metadata': self.metadata,
            'created_at': self.created_at
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'PatientData':
        """Desérialise les donnees patient."""
        patient = cls(data['patient_id'])
        patient.timepoints = data.get('timepoints', {})
        patient.metadata = data.get('metadata', {})
        patient.created_at = data.get('created_at', datetime.now().isoformat())
        return patient


class DataTransformer:
    """
    Transformations de donnees de cytometrie (Logicle, Arcsinh, etc.).
    """
    
    @staticmethod
    def arcsinh_transform(data: np.ndarray, cofactor: float = 150.0) -> np.ndarray:
        """
        Transformation Arcsinh (inverse hyperbolic sine).
        Cofactor standard: 150 pour les donnees CyTOF, 5 pour flow cytometry.
        """
        return np.arcsinh(data / cofactor)
    
    @staticmethod
    def arcsinh_inverse(data: np.ndarray, cofactor: float = 150.0) -> np.ndarray:
        """Inverse de la transformation Arcsinh."""
        return np.sinh(data) * cofactor
    
    @staticmethod
    def logicle_transform(data: np.ndarray, T: float = 262144.0, M: float = 4.5,
                          W: float = 0.5, A: float = 0.0) -> np.ndarray:
        """
        Transformation Logicle (biexponentielle).
        Parametres standard pour les donnees de cytometrie.
        T: Maximum de l'echelle lineaire
        M: Decades de largeur
        W: Linearisation pres de zero
        A: Decades additionnelles (negatifs)
        
        Approximation simplifiee si FlowKit n'est pas disponible.
        """
        if FLOWKIT_AVAILABLE:
            try:
                xform = fk.transforms.LogicleTransform(
                    'logicle', param_t=T, param_m=M, param_w=W, param_a=A
                )
                # FlowKit attend un array 1D par canal
                if len(data.shape) == 1:
                    return xform.apply(data)
                else:
                    result = np.zeros_like(data, dtype=float)
                    for i in range(data.shape[1]):
                        result[:, i] = xform.apply(data[:, i])
                    return result
            except Exception:
                pass
        
        # Approximation si FlowKit absent: Arcsinh modifie
        # Cette approximation capture le comportement log pour les grandes valeurs
        # et lineaire pres de zero
        w_val = W * np.log10(np.e)
        return np.arcsinh(data / (T / (10 ** M))) * (M / np.log(10))
    
    @staticmethod
    def log_transform(data: np.ndarray, base: float = 10.0,
                      min_val: float = 1.0) -> np.ndarray:
        """Transformation logarithmique standard."""
        data_clipped = np.maximum(data, min_val)
        return np.log(data_clipped) / np.log(base)
    
    @staticmethod
    def zscore_normalize(data: np.ndarray) -> np.ndarray:
        """Normalisation Z-score."""
        mean = np.nanmean(data, axis=0)
        std = np.nanstd(data, axis=0)
        std[std == 0] = 1
        return (data - mean) / std
    
    @staticmethod
    def min_max_normalize(data: np.ndarray) -> np.ndarray:
        """Normalisation Min-Max [0, 1]."""
        min_val = np.nanmin(data, axis=0)
        max_val = np.nanmax(data, axis=0)
        range_val = max_val - min_val
        range_val[range_val == 0] = 1
        return (data - min_val) / range_val


class PreGating:
    """
    Pre-gating automatique pour la selection des populations d'interet.
    Basé sur FSC/SSC pour exclure les debris et les doublets.
    Methode robuste avec gestion des valeurs aberrantes et NaN.
    """
    
    @staticmethod
    def find_marker_index(var_names: List[str], patterns: List[str]) -> Optional[int]:
        """Trouve l'index d'un marqueur parmi les patterns donnes."""
        var_upper = [v.upper() for v in var_names]
        for pattern in patterns:
            pattern_upper = pattern.upper()
            for i, name in enumerate(var_upper):
                if pattern_upper in name or name == pattern_upper:
                    return i
        return None
    
    @staticmethod
    def gate_viable_cells(X: np.ndarray, var_names: List[str],
                          min_percentile: float = 2.0, max_percentile: float = 98.0) -> np.ndarray:
        """
        Gate les cellules viables basé sur FSC/SSC.
        Utilise les percentiles pour une detection robuste des outliers.
        
        Args:
            X: Matrice des donnees (n_cells, n_markers)
            var_names: Liste des noms de marqueurs
            min_percentile: Percentile minimum pour exclusion (debris)
            max_percentile: Percentile maximum pour exclusion (doublets)
        
        Returns:
            Masque booléen des cellules viables.
        """
        n_cells = X.shape[0]
        mask = np.ones(n_cells, dtype=bool)
        
        # Trouver FSC (priorite a FSC-A)
        fsc_idx = PreGating.find_marker_index(var_names, ['FSC-A', 'FSC-H', 'FSC'])
        if fsc_idx is not None:
            fsc_vals = X[:, fsc_idx].astype(np.float64)
            # Remplacer les valeurs non-finies par NaN
            fsc_vals = np.where(np.isfinite(fsc_vals), fsc_vals, np.nan)
            
            # Utiliser percentiles pour robustesse
            fsc_min = np.nanpercentile(fsc_vals, min_percentile)
            fsc_max = np.nanpercentile(fsc_vals, max_percentile)
            
            fsc_mask = (fsc_vals >= fsc_min) & (fsc_vals <= fsc_max)
            mask &= np.where(np.isnan(fsc_vals), False, fsc_mask)
        
        # Trouver SSC (priorite a SSC-A)
        ssc_idx = PreGating.find_marker_index(var_names, ['SSC-A', 'SSC-H', 'SSC'])
        if ssc_idx is not None:
            ssc_vals = X[:, ssc_idx].astype(np.float64)
            ssc_vals = np.where(np.isfinite(ssc_vals), ssc_vals, np.nan)
            
            ssc_min = np.nanpercentile(ssc_vals, min_percentile)
            ssc_max = np.nanpercentile(ssc_vals, max_percentile)
            
            ssc_mask = (ssc_vals >= ssc_min) & (ssc_vals <= ssc_max)
            mask &= np.where(np.isnan(ssc_vals), False, ssc_mask)
        
        return mask
    
    @staticmethod
    def gate_singlets(X: np.ndarray, var_names: List[str],
                      ratio_min: float = 0.6, ratio_max: float = 1.5) -> np.ndarray:
        """
        Gate les singlets basé sur le ratio FSC-A/FSC-H.
        Les doublets ont typiquement un ratio > 1.3-1.5.
        
        Args:
            X: Matrice des donnees
            var_names: Liste des noms de marqueurs
            ratio_min: Ratio minimum acceptable (default 0.6)
            ratio_max: Ratio maximum acceptable (default 1.5)
        
        Returns:
            Masque booléen des singlets.
        """
        n_cells = X.shape[0]
        
        fsc_a_idx = PreGating.find_marker_index(var_names, ['FSC-A'])
        fsc_h_idx = PreGating.find_marker_index(var_names, ['FSC-H'])
        
        if fsc_a_idx is None or fsc_h_idx is None:
            # Si pas les deux canaux, retourner True pour toutes les cellules
            return np.ones(n_cells, dtype=bool)
        
        fsc_a = X[:, fsc_a_idx].astype(np.float64)
        fsc_h = X[:, fsc_h_idx].astype(np.float64)
        
        # Gerer les valeurs non-finies
        fsc_a = np.where(np.isfinite(fsc_a), fsc_a, np.nan)
        fsc_h = np.where(np.isfinite(fsc_h), fsc_h, np.nan)
        
        # Eviter division par zero - valeurs tres petites traitees comme invalides
        min_val = 100  # Seuil minimum pour FSC-H valide
        valid_h = fsc_h > min_val
        
        ratio = np.full(n_cells, np.nan)
        ratio[valid_h] = fsc_a[valid_h] / fsc_h[valid_h]
        
        # Masque: ratio valide et dans les limites
        mask = np.isfinite(ratio) & (ratio >= ratio_min) & (ratio <= ratio_max)
        
        return mask
    
    @staticmethod
    def gate_cd45_positive(X: np.ndarray, var_names: List[str],
                           threshold_percentile: float = 10) -> np.ndarray:
        """
        Gate les cellules CD45+ (leucocytes).
        
        Returns:
            Masque booléen des cellules CD45+.
        """
        n_cells = X.shape[0]
        
        cd45_idx = PreGating.find_marker_index(var_names, ['CD45', 'cd45', 'CD45-PECY5', 'CD45-PC5'])
        if cd45_idx is None:
            return np.ones(n_cells, dtype=bool)
        
        cd45_vals = X[:, cd45_idx].astype(np.float64)
        cd45_vals = np.where(np.isfinite(cd45_vals), cd45_vals, np.nan)
        
        threshold = np.nanpercentile(cd45_vals, threshold_percentile)
        
        return np.where(np.isnan(cd45_vals), False, cd45_vals > threshold)


# =============================================================================
# CLASSES CONFORMITÉ ELN - NBM REFERENCE & MRD DETECTION [2][3]
# =============================================================================

class NBMReferenceBuilder:
    """
    Construction d'une reference NBM (Normal Bone Marrow) "frozen" MST.
    Selon [2]: "reference FlowSOM display was built from 19 healthy NBM samples"
    
    Cette classe permet de construire et sauvegarder une reference NBM
    pour comparer les echantillons patients.
    """
    
    def __init__(self, n_nodes: int = 100, seed: int = 42):
        self.n_nodes = n_nodes
        self.seed = seed
        self.frozen_som = None
        self.nbm_node_frequencies: Dict[int, float] = {}
        self.nbm_node_phenotypes: Dict[int, np.ndarray] = {}
        self.nbm_samples_count = 0
        self.is_built = False
        
    def build_reference(self, nbm_data_list: List[np.ndarray], 
                        var_names: List[str]) -> bool:
        """
        Construit la reference FlowSOM frozen depuis plusieurs NBM sains.
        
        Args:
            nbm_data_list: Liste des matrices de donnees NBM
            var_names: Noms des colonnes/marqueurs
            
        Returns:
            True si construction reussie
        """
        if not FLOWSOM_AVAILABLE:
            return False
        
        if len(nbm_data_list) < 15:
            print(f"ATTENTION: {len(nbm_data_list)} NBM < 15 recommandes")
        
        try:
            # Merge tous les NBM
            merged_nbm = np.vstack(nbm_data_list)
            self.nbm_samples_count = len(nbm_data_list)
            
            # Creer AnnData
            adata = ad.AnnData(merged_nbm)
            adata.var_names = var_names
            
            # Train FlowSOM avec n_nodes (standard ISAC: 100)
            xdim = int(np.sqrt(self.n_nodes))
            ydim = self.n_nodes // xdim
            
            self.frozen_som = fs.FlowSOM(
                adata,
                cols_to_use=list(range(len(var_names))),
                xdim=xdim,
                ydim=ydim,
                n_clusters=10,  # Metaclusters
                seed=self.seed
            )
            
            # Calculer les frequences par node dans NBM
            node_assignments = self.frozen_som.get_cell_data().obs['clustering'].values
            total_cells = len(node_assignments)
            
            for node_id in range(self.n_nodes):
                count = np.sum(node_assignments == str(node_id))
                self.nbm_node_frequencies[node_id] = count / total_cells
            
            # Calculer les phenotypes moyens par node
            cell_data = self.frozen_som.get_cell_data()
            X = cell_data.X
            for node_id in range(self.n_nodes):
                mask = node_assignments == str(node_id)
                if mask.sum() > 0:
                    self.nbm_node_phenotypes[node_id] = np.mean(X[mask], axis=0)
                else:
                    self.nbm_node_phenotypes[node_id] = np.zeros(X.shape[1])
            
            self.is_built = True
            return True
            
        except Exception as e:
            print(f"Erreur construction NBM reference: {e}")
            return False
    
    def save_reference(self, filepath: str) -> bool:
        """Sauvegarde la reference NBM sur disque."""
        if not self.is_built:
            return False
        
        try:
            data = {
                'n_nodes': self.n_nodes,
                'seed': self.seed,
                'nbm_samples_count': self.nbm_samples_count,
                'nbm_node_frequencies': self.nbm_node_frequencies,
                'nbm_node_phenotypes': {k: v.tolist() for k, v in self.nbm_node_phenotypes.items()}
            }
            with open(filepath, 'w') as f:
                json.dump(data, f, indent=2)
            return True
        except Exception as e:
            print(f"Erreur sauvegarde NBM reference: {e}")
            return False
    
    def load_reference(self, filepath: str) -> bool:
        """Charge une reference NBM depuis le disque."""
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
            
            self.n_nodes = data['n_nodes']
            self.seed = data['seed']
            self.nbm_samples_count = data['nbm_samples_count']
            self.nbm_node_frequencies = {int(k): v for k, v in data['nbm_node_frequencies'].items()}
            self.nbm_node_phenotypes = {int(k): np.array(v) for k, v in data['nbm_node_phenotypes'].items()}
            self.is_built = True
            return True
        except Exception as e:
            print(f"Erreur chargement NBM reference: {e}")
            return False


class ClonePersistenceTracker:
    """
    Suivi de la persistance clonale T0 -> FU.
    Tracking clone T0 → FU1 → FU2
    """
    
    def __init__(self, nbm_reference: Optional[NBMReferenceBuilder] = None):
        self.nbm_reference = nbm_reference
        self.timepoints: Dict[str, Dict[int, float]] = {}  # {timepoint: {node: freq}}
    
    def add_timepoint(self, timepoint_id: str, node_frequencies: Dict[int, float]):
        """Ajoute un point de temps avec les frequences par node."""
        self.timepoints[timepoint_id] = node_frequencies
    
    def track_clone_persistence(self, diagnosis_id: str, followup_id: str) -> Dict[str, Any]:
        """
        Compare les timepoints pour detecter la persistance clonale.
        
        Args:
            diagnosis_id: ID du timepoint diagnostic (T0)
            followup_id: ID du timepoint de suivi (FU)
            
        Returns:
            Dict avec resultats du tracking
        """
        if diagnosis_id not in self.timepoints or followup_id not in self.timepoints:
            return {'error': 'Timepoint not found'}
        
        diag_freqs = self.timepoints[diagnosis_id]
        fu_freqs = self.timepoints[followup_id]
        
        result = {
            'diagnosis_id': diagnosis_id,
            'followup_id': followup_id,
            'persistent_nodes': [],
            'resolved_nodes': [],
            'new_nodes': []
        }
        
        all_nodes = set(diag_freqs.keys()) | set(fu_freqs.keys())
        
        for node_id in all_nodes:
            diag_freq = diag_freqs.get(node_id, 0)
            fu_freq = fu_freqs.get(node_id, 0)
            
            # Obtenir freq NBM si disponible
            nbm_freq = 0.0
            if self.nbm_reference and self.nbm_reference.is_built:
                nbm_freq = self.nbm_reference.nbm_node_frequencies.get(node_id, 0)
            
            node_data = {
                'node_id': node_id,
                'diagnosis_freq': diag_freq,
                'followup_freq': fu_freq,
                'nbm_freq': nbm_freq,
                'fold_change_to_diag': fu_freq / (diag_freq + 1e-10),
                'fold_change_to_nbm': fu_freq / (nbm_freq + 1e-10) if nbm_freq > 0 else float('inf')
            }
            
            if diag_freq > 0.001 and fu_freq > 0.001:
                # Clone present aux deux timepoints
                node_data['status'] = 'Persistent'
                result['persistent_nodes'].append(node_data)
            elif diag_freq > 0.001 and fu_freq <= 0.001:
                # Clone resolu
                node_data['status'] = 'Resolved'
                result['resolved_nodes'].append(node_data)
            elif diag_freq <= 0.001 and fu_freq > 0.001:
                # Nouveau clone (attention: possible sous-clone emergent)
                node_data['status'] = 'New/Emerging'
                result['new_nodes'].append(node_data)
        
        return result


class PDFReportGenerator:
    """
    Generateur de rapports PDF pour les resultats FlowSOM.
    """
    
    def __init__(self, output_path: str):
        self.output_path = output_path
        self.styles = getSampleStyleSheet() if REPORTLAB_AVAILABLE else None
        
    def generate_report(self, result: Dict[str, Any], figures: List[str],
                        patient_info: Optional[Dict[str, Any]] = None,
                        lsc_result: Optional[Dict[str, Any]] = None) -> bool:
        """
        Genere un rapport PDF complet.
        
        Args:
            result: Resultats FlowSOM
            figures: Liste des chemins vers les figures a inclure
            patient_info: Informations patient optionnelles
            lsc_result: Resultats du score LSC optionnels
            
        Returns:
            True si le rapport a ete genere avec succes.
        """
        if not REPORTLAB_AVAILABLE:
            return False
        
        try:
            doc = SimpleDocTemplate(
                self.output_path,
                pagesize=A4,
                rightMargin=2*cm,
                leftMargin=2*cm,
                topMargin=2*cm,
                bottomMargin=2*cm
            )
            
            story = []
            
            # Titre
            title_style = ParagraphStyle(
                'CustomTitle',
                parent=self.styles['Heading1'],
                fontSize=24,
                textColor=HexColor('#1e3a5f'),
                spaceAfter=20,
                alignment=1  # Centré
            )
            story.append(Paragraph("Rapport d'Analyse FlowSOM", title_style))
            story.append(Spacer(1, 0.5*cm))
            
            # Date et heure
            date_style = ParagraphStyle(
                'DateStyle',
                parent=self.styles['Normal'],
                fontSize=10,
                textColor=HexColor('#666666'),
                alignment=1
            )
            story.append(Paragraph(
                f"Genere le {datetime.now().strftime('%d/%m/%Y a %H:%M')}",
                date_style
            ))
            story.append(Spacer(1, 1*cm))
            
            # Informations patient si disponibles
            if patient_info:
                story.append(Paragraph("Informations Patient", self.styles['Heading2']))
                patient_data = [
                    ['ID Patient:', patient_info.get('patient_id', 'N/A')],
                    ['Date prelèvement:', patient_info.get('date', 'N/A')],
                    ['Notes:', patient_info.get('notes', '')]
                ]
                table = Table(patient_data, colWidths=[4*cm, 12*cm])
                table.setStyle(TableStyle([
                    ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
                    ('FONTSIZE', (0, 0), (-1, -1), 10),
                    ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
                    ('ALIGN', (0, 0), (0, -1), 'RIGHT'),
                    ('ALIGN', (1, 0), (1, -1), 'LEFT'),
                    ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
                ]))
                story.append(table)
                story.append(Spacer(1, 0.5*cm))
            
            # Parametres d'analyse
            story.append(Paragraph("Parametres d'Analyse", self.styles['Heading2']))
            params = result.get('params', {})
            params_data = [
                ['Grille SOM:', f"{params.get('xdim', 10)} x {params.get('ydim', 10)}"],
                ['Metaclusters:', str(result.get('n_clusters', 'N/A'))],
                ['Nombre de cellules:', f"{result.get('cell_data').shape[0]:,}" if result.get('cell_data') is not None else 'N/A'],
                ['Marqueurs utilises:', str(len(result.get('cols_to_use', [])))]
            ]
            table = Table(params_data, colWidths=[5*cm, 11*cm])
            table.setStyle(TableStyle([
                ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 10),
                ('BACKGROUND', (0, 0), (-1, -1), HexColor('#f5f5f5')),
                ('BOX', (0, 0), (-1, -1), 1, colors.lightgrey),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
                ('TOPPADDING', (0, 0), (-1, -1), 8),
            ]))
            story.append(table)
            story.append(Spacer(1, 0.5*cm))
            
            # Score LSC si disponible
            if lsc_result:
                story.append(Paragraph("Score LSC (Leukemic Stem Cell)", self.styles['Heading2']))
                lsc_data = [
                    ['Pourcentage LSC:', f"{lsc_result.get('lsc_percentage', 0):.2f}%"],
                    ['CD34+:', f"{lsc_result.get('cd34_pos_percentage', 0):.2f}%"],
                    ['CD38- (parmi CD34+):', f"{lsc_result.get('cd38_neg_percentage', 0):.2f}%"],
                    ['Marqueurs detectes:', ', '.join(lsc_result.get('available_markers', []))],
                    ['Marqueurs manquants:', ', '.join(lsc_result.get('missing_markers', []))]
                ]
                table = Table(lsc_data, colWidths=[5*cm, 11*cm])
                table.setStyle(TableStyle([
                    ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
                    ('FONTSIZE', (0, 0), (-1, -1), 10),
                    ('BACKGROUND', (0, 0), (-1, 0), HexColor('#e3f2fd')),
                    ('BOX', (0, 0), (-1, -1), 1, colors.lightgrey),
                    ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
                    ('TOPPADDING', (0, 0), (-1, -1), 8),
                ]))
                story.append(table)
                story.append(Spacer(1, 0.5*cm))
            
            # Figures
            if figures:
                story.append(PageBreak())
                story.append(Paragraph("Visualisations", self.styles['Heading2']))
                story.append(Spacer(1, 0.5*cm))
                
                for fig_path in figures:
                    if os.path.exists(fig_path):
                        img = Image(fig_path, width=16*cm, height=12*cm)
                        story.append(img)
                        story.append(Spacer(1, 0.5*cm))
            
            doc.build(story)
            return True
            
        except Exception as e:
            print(f"Erreur generation PDF: {e}")
            return False


# =============================================================================
# WORKER THREAD POUR LE CALCUL FLOWSOM
# =============================================================================

class FlowSOMWorker(QThread):
    """
    Thread de travail pour exécuter l'analyse FlowSOM sans bloquer l'interface.
    """
    
    progress = pyqtSignal(int)
    status = pyqtSignal(str)
    finished = pyqtSignal(object)
    error = pyqtSignal(str)
    
    def __init__(
        self,
        healthy_files: List[str],
        pathological_files: List[str],
        params: Dict[str, Any]
    ):
        super().__init__()
        self.healthy_files = healthy_files
        self.pathological_files = pathological_files
        self.params = params
        self._is_cancelled = False
        
    def cancel(self):
        self._is_cancelled = True
        
    def run(self):
        try:
            # Etape 1: Chargement des fichiers sains
            self.status.emit("Chargement des fichiers sains...")
            self.progress.emit(5)
            
            healthy_data_list = []
            for fpath in self.healthy_files:
                if self._is_cancelled:
                    return
                try:
                    adata = fs.io.read_FCS(fpath)
                    adata.obs['condition'] = 'Sain'
                    adata.obs['file_origin'] = Path(fpath).name
                    healthy_data_list.append(adata)
                except Exception as e:
                    self.status.emit(f"Erreur: {Path(fpath).name}")
                    
            self.progress.emit(20)
            
            # Etape 2: Chargement des fichiers pathologiques
            self.status.emit("Chargement des fichiers pathologiques...")
            
            pathological_data_list = []
            for fpath in self.pathological_files:
                if self._is_cancelled:
                    return
                try:
                    adata = fs.io.read_FCS(fpath)
                    adata.obs['condition'] = 'Pathologique'
                    adata.obs['file_origin'] = Path(fpath).name
                    pathological_data_list.append(adata)
                except Exception as e:
                    self.status.emit(f"Erreur: {Path(fpath).name}")
                    
            self.progress.emit(40)
            
            if not healthy_data_list and not pathological_data_list:
                self.error.emit("Aucun fichier FCS valide n'a pu etre charge.")
                return
                
            # Etape 3: Concatenation
            self.status.emit("Concatenation des donnees...")
            
            all_data = healthy_data_list + pathological_data_list
            
            if len(all_data) > 1:
                common_vars = set(all_data[0].var_names)
                for adata in all_data[1:]:
                    common_vars &= set(adata.var_names)
                common_vars = list(common_vars)
                all_data = [adata[:, common_vars].copy() for adata in all_data]
            
            combined_data = ad.concat(all_data, join='inner')
            self.progress.emit(50)
            
            # Etape 4: Pre-gating (optionnel)
            if self.params.get('pregate', False):
                self.status.emit("Application du pre-gating...")
                X = combined_data.X
                if hasattr(X, 'toarray'):
                    X = X.toarray()
                
                var_names_list = list(combined_data.var_names)
                
                # Gate cellules viables
                viable_mask = PreGating.gate_viable_cells(X, var_names_list)
                # Gate singlets
                singlet_mask = PreGating.gate_singlets(X, var_names_list)
                # Combiner
                final_mask = viable_mask & singlet_mask
                
                n_before = combined_data.shape[0]
                combined_data = combined_data[final_mask, :].copy()
                n_after = combined_data.shape[0]
                
                self.status.emit(f"Pre-gating: {n_before} -> {n_after} cellules ({n_after/n_before*100:.1f}%)")
            
            # Etape 5: Transformation des donnees (optionnel)
            transform_type = self.params.get('transform', 'Aucune')
            if transform_type != 'Aucune':
                self.status.emit(f"Application transformation {transform_type}...")
                X = combined_data.X
                if hasattr(X, 'toarray'):
                    X = X.toarray()
                
                cofactor = self.params.get('cofactor', 5.0)
                
                if 'Arcsinh' in transform_type:
                    if 'cofactor=5' in transform_type:
                        cofactor = 5.0
                    elif 'cofactor=150' in transform_type:
                        cofactor = 150.0
                    X_transformed = DataTransformer.arcsinh_transform(X, cofactor)
                elif transform_type == 'Logicle':
                    X_transformed = DataTransformer.logicle_transform(X)
                elif transform_type == 'Log10':
                    X_transformed = DataTransformer.log_transform(X)
                elif transform_type == 'Z-score':
                    X_transformed = DataTransformer.zscore_normalize(X)
                else:
                    X_transformed = X
                
                combined_data.X = X_transformed
            
            # Etape 6: Colonnes a utiliser
            self.status.emit("Preparation des parametres...")
            
            var_names = list(combined_data.var_names)
            
            if self.params.get('exclude_scatter', True):
                exclude_patterns = ['Time', 'FSC', 'SSC', 'Event', 'event']
                cols_to_use = []
                for i, name in enumerate(var_names):
                    if not any(pattern.lower() in name.lower() for pattern in exclude_patterns):
                        cols_to_use.append(i)
            else:
                cols_to_use = list(range(len(var_names)))
            
            if len(cols_to_use) == 0:
                cols_to_use = list(range(len(var_names)))
                
            self.progress.emit(55)
            
            # Etape 5: FlowSOM
            self.status.emit("Execution de FlowSOM...")
            
            xdim = self.params.get('xdim', 10)
            ydim = self.params.get('ydim', 10)
            n_clusters = self.params.get('n_clusters', 10)
            seed = self.params.get('seed', 42)
            
            if self.params.get('auto_cluster', False):
                self.status.emit("Recherche du nombre optimal de clusters...")
                n_clusters = self._find_optimal_clusters(combined_data, cols_to_use, xdim, ydim, seed)
                self.status.emit(f"Nombre optimal: {n_clusters} clusters")
            
            self.progress.emit(60)
            
            fsom = fs.FlowSOM(
                combined_data,
                cols_to_use=cols_to_use,
                xdim=xdim,
                ydim=ydim,
                n_clusters=n_clusters,
                seed=seed
            )
            
            self.progress.emit(90)
            
            # Etape 6: Resultats
            self.status.emit("Preparation des resultats...")
            
            cell_data = fsom.get_cell_data()
            cell_data.obs['condition'] = combined_data.obs['condition'].values
            cell_data.obs['file_origin'] = combined_data.obs['file_origin'].values
            
            result = {
                'fsom': fsom,
                'combined_data': combined_data,
                'cell_data': cell_data,
                'cluster_data': fsom.get_cluster_data(),
                'n_clusters': n_clusters,
                'cols_to_use': cols_to_use,
                'var_names': var_names,
                'params': self.params
            }
            
            self.progress.emit(100)
            self.status.emit("Analyse terminee avec succes!")
            self.finished.emit(result)
            
        except Exception as e:
            import traceback
            self.error.emit(f"Erreur: {str(e)}\n{traceback.format_exc()}")
            
    def _find_optimal_clusters(self, data, cols_to_use, xdim, ydim, seed, max_clusters=20):
        n_samples = min(5000, data.shape[0])
        idx = np.random.choice(data.shape[0], n_samples, replace=False)
        subset_data = data[idx, :].copy()
        
        X = subset_data.X[:, cols_to_use]
        if hasattr(X, 'toarray'):
            X = X.toarray()
        X = np.nan_to_num(X, nan=0.0)
        
        scores = []
        cluster_range = range(2, min(max_clusters + 1, n_samples // 10))
        
        for k in cluster_range:
            try:
                clustering = AgglomerativeClustering(n_clusters=k)
                labels = clustering.fit_predict(X)
                score = silhouette_score(X, labels, sample_size=min(1000, n_samples))
                scores.append((k, score))
            except Exception:
                continue
                
        return max(scores, key=lambda x: x[1])[0] if scores else 10


# =============================================================================
# WORKER THREAD POUR UMAP / t-SNE (evite le freeze de l'interface)
# =============================================================================

class EmbeddingWorker(QThread):
    """
    Thread de travail pour calculer UMAP ou t-SNE sans bloquer l'interface.
    """
    
    progress = pyqtSignal(int)
    status = pyqtSignal(str)
    finished = pyqtSignal(object)
    error = pyqtSignal(str)
    
    def __init__(self, cell_data, cols_to_use: List[int], method: str = 'umap', 
                 max_cells: int = 5000, n_neighbors: int = 15, min_dist: float = 0.1,
                 perplexity: int = 30):
        super().__init__()
        self.cell_data = cell_data
        self.cols_to_use = cols_to_use
        self.method = method  # 'umap' ou 'tsne'
        self.max_cells = max_cells
        self.n_neighbors = n_neighbors
        self.min_dist = min_dist
        self.perplexity = perplexity
        self._is_cancelled = False
        
    def cancel(self):
        self._is_cancelled = True
        
    def run(self):
        try:
            self.status.emit(f"Preparation des donnees pour {self.method.upper()}...")
            self.progress.emit(10)
            
            # Extraire les donnees
            X = self.cell_data.X
            if hasattr(X, 'toarray'):
                X = X.toarray()
            
            # Selectionner les colonnes utilisees
            X_subset = X[:, self.cols_to_use]
            
            # Nettoyer les NaN/Inf
            X_subset = np.nan_to_num(X_subset, nan=0.0, posinf=0.0, neginf=0.0)
            
            # Sous-echantillonnage pour la performance
            n_cells = X_subset.shape[0]
            if n_cells > self.max_cells:
                self.status.emit(f"Sous-echantillonnage: {self.max_cells}/{n_cells} cellules...")
                idx = np.random.choice(n_cells, self.max_cells, replace=False)
                X_sample = X_subset[idx, :]
                metaclusters_sample = self.cell_data.obs['metaclustering'].values[idx]
                conditions_sample = self.cell_data.obs['condition'].values[idx]
            else:
                X_sample = X_subset
                idx = np.arange(n_cells)
                metaclusters_sample = self.cell_data.obs['metaclustering'].values
                conditions_sample = self.cell_data.obs['condition'].values
            
            self.progress.emit(30)
            
            if self._is_cancelled:
                return
            
            # Calcul de l'embedding
            if self.method == 'umap':
                self.status.emit("Calcul UMAP en cours (peut prendre 30-60s)...")
                self.progress.emit(50)
                
                if UMAP_AVAILABLE:
                    reducer = umap.UMAP(
                        n_neighbors=self.n_neighbors,
                        min_dist=self.min_dist,
                        n_components=2,
                        metric='euclidean',
                        random_state=42
                    )
                    embedding = reducer.fit_transform(X_sample)
                elif SCANPY_AVAILABLE:
                    # Fallback sur scanpy
                    adata_temp = ad.AnnData(X_sample)
                    sc.pp.neighbors(adata_temp, n_neighbors=self.n_neighbors)
                    sc.tl.umap(adata_temp, min_dist=self.min_dist)
                    embedding = adata_temp.obsm['X_umap']
                else:
                    raise RuntimeError("Ni umap-learn ni scanpy ne sont installes")
                    
            else:  # t-SNE
                self.status.emit("Calcul t-SNE en cours (peut prendre 1-2min)...")
                self.progress.emit(50)
                
                if TSNE_AVAILABLE:
                    # Ajuster perplexity si necessaire
                    perp = min(self.perplexity, X_sample.shape[0] // 4)
                    perp = max(5, perp)
                    
                    tsne = TSNE(
                        n_components=2,
                        perplexity=perp,
                        learning_rate='auto',
                        init='pca',
                        random_state=42,
                        n_iter=1000
                    )
                    embedding = tsne.fit_transform(X_sample)
                elif SCANPY_AVAILABLE:
                    adata_temp = ad.AnnData(X_sample)
                    sc.tl.tsne(adata_temp, perplexity=min(self.perplexity, X_sample.shape[0] // 4))
                    embedding = adata_temp.obsm['X_tsne']
                else:
                    raise RuntimeError("Ni sklearn ni scanpy ne sont installes")
            
            if self._is_cancelled:
                return
                
            self.progress.emit(90)
            self.status.emit("Finalisation...")
            
            result = {
                'embedding': embedding,
                'metaclusters': metaclusters_sample,
                'conditions': conditions_sample,
                'indices': idx,
                'method': self.method
            }
            
            self.progress.emit(100)
            self.status.emit(f"{self.method.upper()} termine!")
            self.finished.emit(result)
            
        except Exception as e:
            import traceback
            self.error.emit(f"Erreur {self.method.upper()}: {str(e)}\n{traceback.format_exc()}")


# =============================================================================
# WIDGET DE VISUALISATION MATPLOTLIB
# =============================================================================

class MatplotlibCanvas(FigureCanvas):
    def __init__(self, parent=None, width=8, height=6, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.patch.set_facecolor('#1e1e2e')
        super().__init__(self.fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.updateGeometry()
        
    def clear_figure(self):
        self.fig.clear()
        self.fig.patch.set_facecolor('#1e1e2e')
        self.draw()


# =============================================================================
# FENÊTRE PRINCIPALE DE L'APPLICATION
# =============================================================================

class FlowSOMApp(QMainWindow):
    def __init__(self):
        super().__init__()
        
        # Configuration chargee depuis config.yaml
        self.config = CONFIG
        
        self.healthy_folder: Optional[str] = None
        self.pathological_folder: Optional[str] = None
        self.healthy_files: List[str] = []
        self.pathological_files: List[str] = []
        self.result: Optional[Dict[str, Any]] = None
        self.worker: Optional[FlowSOMWorker] = None
        self.embedding_worker: Optional[EmbeddingWorker] = None
        self.embedding_result: Optional[Dict[str, Any]] = None
        
        # Attributs patient
        self.current_patient: Optional[PatientData] = None
        self.patients_db: Dict[str, PatientData] = {}
        self.transformed_data: Optional[np.ndarray] = None
        
        # Reference NBM
        self.nbm_reference: Optional[NBMReferenceBuilder] = None
        self.clone_tracker: Optional[ClonePersistenceTracker] = None
        self.nbm_files: List[str] = []  # Fichiers NBM pour reference
        
        # Fichier FCS charge pour visualisation
        self.current_fcs_adata = None
        
        self.setWindowTitle("FlowSOM Analyzer - Analyse Cytométrie")
        self.setMinimumSize(1500, 950)
        self.setStyleSheet(STYLESHEET)
        
        self._create_ui()
        self.statusBar().showMessage("Pret. Chargez des fichiers FCS pour commencer l'analyse.")
        
    def _create_ui(self):
        central_widget = QWidget()
        central_widget.setObjectName("centralWidget")
        self.setCentralWidget(central_widget)
        
        main_layout = QHBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)
        
        # =====================================================================
        # PANNEAU DE GAUCHE
        # =====================================================================
        left_panel = QWidget()
        left_panel.setObjectName("leftPanel")
        left_panel.setFixedWidth(380)
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(20, 25, 20, 20)
        left_layout.setSpacing(15)
        
        # Header avec titre
        header_widget = QWidget()
        header_layout = QVBoxLayout(header_widget)
        header_layout.setSpacing(5)
        header_layout.setContentsMargins(0, 0, 0, 15)
        
        # Icone et titre
        title_layout = QHBoxLayout()
        icon_label = QLabel("[F]")
        icon_label.setStyleSheet("font-size: 24pt; background: transparent; color: #89b4fa; font-weight: bold;")
        title_layout.addWidget(icon_label)
        
        title_text = QLabel("FlowSOM Analyzer")
        title_text.setStyleSheet("""
            font-size: 20pt;
            font-weight: 700;
            color: #89b4fa;
            background: transparent;
        """)
        title_layout.addWidget(title_text)
        title_layout.addStretch()
        header_layout.addLayout(title_layout)
        
        subtitle = QLabel("Analyse FlowSOM - Reference NBM vs Patient Diagnostic")
        subtitle.setObjectName("subtitleLabel")
        header_layout.addWidget(subtitle)
        
        # Ligne décorative
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setStyleSheet("background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #89b4fa, stop:0.5 #cba6f7, stop:1 transparent); min-height: 2px; border: none;")
        header_layout.addWidget(line)
        
        left_layout.addWidget(header_widget)
        
        # Scroll area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        scroll_layout.setSpacing(18)
        scroll_layout.setContentsMargins(0, 0, 10, 0)
        
        # ----- GROUPE: Chargement des donnees -----
        data_group = QGroupBox("Chargement des Donnees")
        data_layout = QVBoxLayout(data_group)
        data_layout.setSpacing(12)
        
        # Mode d'analyse
        mode_layout = QHBoxLayout()
        lbl_mode = QLabel("Mode:")
        lbl_mode.setObjectName("sectionLabel")
        mode_layout.addWidget(lbl_mode)
        
        self.chk_compare_nbm = QCheckBox("Comparer avec NBM")
        self.chk_compare_nbm.setChecked(False)
        self.chk_compare_nbm.setToolTip("Cochez pour comparer Patient vs Reference NBM. Decochez pour analyser un patient seul.")
        self.chk_compare_nbm.stateChanged.connect(self._toggle_nbm_mode)
        mode_layout.addWidget(self.chk_compare_nbm)
        mode_layout.addStretch()
        data_layout.addLayout(mode_layout)
        
        # Reference NBM (Normal Bone Marrow)
        self.btn_healthy = QPushButton("Reference NBM (Moelle Normale)")
        self.btn_healthy.setObjectName("successBtn")
        self.btn_healthy.setCursor(Qt.PointingHandCursor)
        self.btn_healthy.clicked.connect(self._load_healthy_folder)
        self.btn_healthy.setToolTip("Fichier ou dossier de reference NBM sain")
        self.btn_healthy.setEnabled(False)  # Desactive par defaut
        data_layout.addWidget(self.btn_healthy)
        
        self.lbl_healthy = QLabel("Aucune reference selectionnee")
        self.lbl_healthy.setObjectName("fileLabel")
        self.lbl_healthy.setWordWrap(True)
        data_layout.addWidget(self.lbl_healthy)
        
        # Patient Diagnostic
        self.btn_pathological = QPushButton("Patient Diagnostic")
        self.btn_pathological.setObjectName("dangerBtn")
        self.btn_pathological.setCursor(Qt.PointingHandCursor)
        self.btn_pathological.clicked.connect(self._load_pathological_folder)
        self.btn_pathological.setToolTip("Fichier ou dossier patient a analyser")
        data_layout.addWidget(self.btn_pathological)
        
        self.lbl_pathological = QLabel("Aucun patient selectionne")
        self.lbl_pathological.setObjectName("fileLabel")
        self.lbl_pathological.setWordWrap(True)
        data_layout.addWidget(self.lbl_pathological)
        
        # Info fichiers
        self.lbl_file_info = QLabel("Fichiers: 0 reference, 0 patient")
        self.lbl_file_info.setStyleSheet("color: #a6adc8; font-weight: 500; padding: 5px;")
        data_layout.addWidget(self.lbl_file_info)
        
        scroll_layout.addWidget(data_group)
        
        # ----- GROUPE: Parametres FlowSOM -----
        params_group = QGroupBox("Parametres FlowSOM")
        params_layout = QGridLayout(params_group)
        params_layout.setSpacing(12)
        params_layout.setColumnStretch(1, 1)
        
        row = 0
        
        # X dimension
        lbl_xdim = QLabel("X Dimension:")
        lbl_xdim.setObjectName("sectionLabel")
        params_layout.addWidget(lbl_xdim, row, 0)
        self.spin_xdim = QSpinBox()
        self.spin_xdim.setRange(5, 50)
        self.spin_xdim.setValue(10)
        self.spin_xdim.setToolTip("Dimension X de la grille SOM")
        params_layout.addWidget(self.spin_xdim, row, 1)
        row += 1
        
        # Y dimension
        lbl_ydim = QLabel("Y Dimension:")
        lbl_ydim.setObjectName("sectionLabel")
        params_layout.addWidget(lbl_ydim, row, 0)
        self.spin_ydim = QSpinBox()
        self.spin_ydim.setRange(5, 50)
        self.spin_ydim.setValue(10)
        self.spin_ydim.setToolTip("Dimension Y de la grille SOM")
        params_layout.addWidget(self.spin_ydim, row, 1)
        row += 1
        
        # Metaclusters
        lbl_clusters = QLabel("Métaclusters:")
        lbl_clusters.setObjectName("sectionLabel")
        params_layout.addWidget(lbl_clusters, row, 0)
        self.spin_n_clusters = QSpinBox()
        self.spin_n_clusters.setRange(2, 50)
        self.spin_n_clusters.setValue(10)
        self.spin_n_clusters.setToolTip("Nombre de métaclusters")
        params_layout.addWidget(self.spin_n_clusters, row, 1)
        row += 1
        
        # Seed
        lbl_seed = QLabel("Seed:")
        lbl_seed.setObjectName("sectionLabel")
        params_layout.addWidget(lbl_seed, row, 0)
        self.spin_seed = QSpinBox()
        self.spin_seed.setRange(0, 99999)
        self.spin_seed.setValue(42)
        self.spin_seed.setToolTip("Graine pour reproductibilité")
        params_layout.addWidget(self.spin_seed, row, 1)
        row += 1
        
        # Checkboxes
        self.chk_auto_cluster = QCheckBox("Auto-clustering (silhouette)")
        self.chk_auto_cluster.setToolTip("Determiner automatiquement le nombre optimal")
        self.chk_auto_cluster.stateChanged.connect(self._toggle_auto_cluster)
        params_layout.addWidget(self.chk_auto_cluster, row, 0, 1, 2)
        row += 1
        
        self.chk_exclude_scatter = QCheckBox("Exclure FSC/SSC/Time")
        self.chk_exclude_scatter.setChecked(True)
        self.chk_exclude_scatter.setToolTip("Exclure les parametres de scatter")
        params_layout.addWidget(self.chk_exclude_scatter, row, 0, 1, 2)
        row += 1
        
        # Pre-gating automatique
        self.chk_pregate = QCheckBox("Pre-gating auto (debris/doublets)")
        self.chk_pregate.setChecked(False)
        self.chk_pregate.setToolTip("Exclure automatiquement debris et doublets")
        params_layout.addWidget(self.chk_pregate, row, 0, 1, 2)
        
        scroll_layout.addWidget(params_group)
        
        # ----- GROUPE: Transformation ----- (conditionnel selon config)
        if self.config['modules_optionnels'].get('data_transformations', True):
            transform_group = QGroupBox("Transformation des donnees")
            transform_layout = QGridLayout(transform_group)
            transform_layout.setSpacing(12)
            
            lbl_transform = QLabel("Type:")
            lbl_transform.setObjectName("sectionLabel")
            transform_layout.addWidget(lbl_transform, 0, 0)
            
            self.combo_transform = QComboBox()
            self.combo_transform.addItems([
                "Aucune",
                "Arcsinh (cofactor=5)",
                "Arcsinh (cofactor=150)",
                "Logicle",
                "Log10",
                "Z-score"
            ])
            self.combo_transform.setToolTip("Transformation appliquee aux donnees")
            transform_layout.addWidget(self.combo_transform, 0, 1)
            
            lbl_cofactor = QLabel("Cofacteur:")
            lbl_cofactor.setObjectName("sectionLabel")
            transform_layout.addWidget(lbl_cofactor, 1, 0)
            
            self.spin_cofactor = QDoubleSpinBox()
            self.spin_cofactor.setRange(1.0, 500.0)
            self.spin_cofactor.setValue(5.0)
            self.spin_cofactor.setToolTip("Cofacteur pour Arcsinh (5 flow, 150 CyTOF)")
            transform_layout.addWidget(self.spin_cofactor, 1, 1)
            
            scroll_layout.addWidget(transform_group)
        else:
            # Widgets par defaut si module desactive
            self.combo_transform = None
            self.spin_cofactor = None
        
        # ----- GROUPE: Reference NBM (ELN) ----- (conditionnel selon config)
        if self.config['modules_avances'].get('nbm_reference', True):
            nbm_group = QGroupBox("Reference NBM (ELN)")
            nbm_layout = QVBoxLayout(nbm_group)
            nbm_layout.setSpacing(10)
            
            self.btn_load_nbm = QPushButton("Charger dossier NBM (reference)")
            self.btn_load_nbm.setToolTip("15-20 echantillons NBM sains pour reference")
            self.btn_load_nbm.clicked.connect(self._load_nbm_folder)
            nbm_layout.addWidget(self.btn_load_nbm)
            
            self.lbl_nbm_status = QLabel("Reference NBM: Non chargee")
            self.lbl_nbm_status.setStyleSheet("color: #a6adc8; font-size: 9pt; padding: 5px;")
            nbm_layout.addWidget(self.lbl_nbm_status)
            
            self.btn_build_nbm_ref = QPushButton("Construire reference frozen MST")
            self.btn_build_nbm_ref.setEnabled(False)
            self.btn_build_nbm_ref.clicked.connect(self._build_nbm_reference)
            nbm_layout.addWidget(self.btn_build_nbm_ref)
            
            self.btn_load_nbm_ref = QPushButton("Charger reference existante")
            self.btn_load_nbm_ref.clicked.connect(self._load_nbm_reference)
            nbm_layout.addWidget(self.btn_load_nbm_ref)
            
            scroll_layout.addWidget(nbm_group)
        else:
            self.btn_load_nbm = None
            self.lbl_nbm_status = None
            self.btn_build_nbm_ref = None
            self.btn_load_nbm_ref = None
        
        # ----- GROUPE: Patient / Time-Series ----- (conditionnel selon config)
        if self.config['modules_optionnels'].get('patient_tracking', True):
            patient_group = QGroupBox("Suivi Patient (Time-Series)")
            patient_layout = QGridLayout(patient_group)
            patient_layout.setSpacing(10)
            
            lbl_patient_id = QLabel("ID Patient:")
            lbl_patient_id.setObjectName("sectionLabel")
            patient_layout.addWidget(lbl_patient_id, 0, 0)
            
            self.edit_patient_id = QLineEdit()
            self.edit_patient_id.setPlaceholderText("Ex: PAT-001")
            patient_layout.addWidget(self.edit_patient_id, 0, 1)
            
            lbl_timepoint = QLabel("Point temps:")
            lbl_timepoint.setObjectName("sectionLabel")
            patient_layout.addWidget(lbl_timepoint, 1, 0)
            
            self.edit_timepoint = QLineEdit()
            self.edit_timepoint.setPlaceholderText("Ex: T0, J30, M3...")
            patient_layout.addWidget(self.edit_timepoint, 1, 1)
            
            self.btn_save_patient = QPushButton("Sauvegarder timepoint")
            self.btn_save_patient.setEnabled(False)
            self.btn_save_patient.clicked.connect(self._save_patient_timepoint)
            patient_layout.addWidget(self.btn_save_patient, 2, 0, 1, 2)
            
            self.btn_load_patient = QPushButton("Charger historique patient")
            self.btn_load_patient.clicked.connect(self._load_patient_history)
            patient_layout.addWidget(self.btn_load_patient, 3, 0, 1, 2)
            
            # Clone tracking (si active)
            if self.config['modules_avances'].get('clone_persistence', True):
                self.btn_track_clones = QPushButton("Analyser persistance clonale")
                self.btn_track_clones.setEnabled(False)
                self.btn_track_clones.setToolTip("Compare T0 vs FU pour detecter clones persistants [2]")
                self.btn_track_clones.clicked.connect(self._analyze_clone_persistence)
                patient_layout.addWidget(self.btn_track_clones, 4, 0, 1, 2)
            else:
                self.btn_track_clones = None
            
            scroll_layout.addWidget(patient_group)
        else:
            self.edit_patient_id = None
            self.edit_timepoint = None
            self.btn_save_patient = None
            self.btn_load_patient = None
            self.btn_track_clones = None
        
        # ----- GROUPE: Actions -----
        actions_group = QGroupBox("Actions")
        actions_layout = QVBoxLayout(actions_group)
        actions_layout.setSpacing(12)
        
        self.btn_run = QPushButton("Lancer l'Analyse FlowSOM")
        self.btn_run.setObjectName("primaryBtn")
        self.btn_run.setEnabled(False)
        self.btn_run.setCursor(Qt.PointingHandCursor)
        self.btn_run.setMinimumHeight(50)
        self.btn_run.clicked.connect(self._run_analysis)
        actions_layout.addWidget(self.btn_run)
        
        self.btn_cancel = QPushButton("Annuler")
        self.btn_cancel.setEnabled(False)
        self.btn_cancel.setCursor(Qt.PointingHandCursor)
        self.btn_cancel.clicked.connect(self._cancel_analysis)
        actions_layout.addWidget(self.btn_cancel)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setFormat("%p%")
        actions_layout.addWidget(self.progress_bar)
        
        self.lbl_status = QLabel("En attente...")
        self.lbl_status.setWordWrap(True)
        self.lbl_status.setAlignment(Qt.AlignCenter)
        self.lbl_status.setStyleSheet("color: #a6adc8; padding: 8px; font-size: 9pt;")
        actions_layout.addWidget(self.lbl_status)
        
        scroll_layout.addWidget(actions_group)
        
        # ----- GROUPE: Visualisation -----
        viz_group = QGroupBox("Visualisation")
        viz_layout = QVBoxLayout(viz_group)
        viz_layout.setSpacing(12)
        
        lbl_viz = QLabel("Type de vue:")
        lbl_viz.setObjectName("sectionLabel")
        viz_layout.addWidget(lbl_viz)
        
        self.combo_viz_type = QComboBox()
        self.combo_viz_type.addItems([
            "Heatmap",
            "Distribution NBM/Patient",
            "Cluster Numbers",
            "Arbre MST (matplotlib)",
            "Grille SOM (matplotlib)",
            "Grid View",
            "Star Chart (MST)",
            "Marker Expression",
            "UMAP Metaclusters",
            "t-SNE Metaclusters"
        ])
        self.combo_viz_type.currentIndexChanged.connect(self._update_visualization)
        viz_layout.addWidget(self.combo_viz_type)
        
        lbl_marker = QLabel("Marqueur:")
        lbl_marker.setObjectName("sectionLabel")
        viz_layout.addWidget(lbl_marker)
        
        self.combo_marker = QComboBox()
        self.combo_marker.setEnabled(False)
        self.combo_marker.currentIndexChanged.connect(self._update_visualization)
        viz_layout.addWidget(self.combo_marker)
        
        # Nombre de cellules pour UMAP/t-SNE
        lbl_n_cells = QLabel("Cellules UMAP/t-SNE:")
        lbl_n_cells.setObjectName("sectionLabel")
        viz_layout.addWidget(lbl_n_cells)
        
        self.spin_embedding_cells = QSpinBox()
        self.spin_embedding_cells.setRange(500, 100000)
        self.spin_embedding_cells.setValue(5000)
        self.spin_embedding_cells.setSingleStep(1000)
        self.spin_embedding_cells.setToolTip("Nombre max de cellules pour UMAP/t-SNE (plus = plus lent)")
        viz_layout.addWidget(self.spin_embedding_cells)
        
        # Checkbox pour utiliser toutes les cellules
        self.chk_all_cells = QCheckBox("Toutes les cellules (lent)")
        self.chk_all_cells.setChecked(False)
        self.chk_all_cells.setToolTip("Utiliser toutes les cellules pour UMAP/t-SNE (peut etre tres lent)")
        self.chk_all_cells.stateChanged.connect(self._toggle_all_cells)
        viz_layout.addWidget(self.chk_all_cells)
        
        scroll_layout.addWidget(viz_group)
        
        # ----- GROUPE: Export -----
        export_group = QGroupBox("Export")
        export_layout = QVBoxLayout(export_group)
        export_layout.setSpacing(10)
        
        self.btn_export_fcs = QPushButton("Exporter en FCS")
        self.btn_export_fcs.setObjectName("exportBtn")
        self.btn_export_fcs.setEnabled(False)
        self.btn_export_fcs.setCursor(Qt.PointingHandCursor)
        self.btn_export_fcs.clicked.connect(self._export_fcs)
        export_layout.addWidget(self.btn_export_fcs)
        
        self.btn_export_csv = QPushButton("Exporter en CSV")
        self.btn_export_csv.setObjectName("exportBtn")
        self.btn_export_csv.setEnabled(False)
        self.btn_export_csv.setCursor(Qt.PointingHandCursor)
        self.btn_export_csv.clicked.connect(self._export_csv)
        export_layout.addWidget(self.btn_export_csv)
        
        self.btn_export_fig = QPushButton("Exporter la Figure")
        self.btn_export_fig.setObjectName("exportBtn")
        self.btn_export_fig.setEnabled(False)
        self.btn_export_fig.setCursor(Qt.PointingHandCursor)
        self.btn_export_fig.clicked.connect(self._export_figure)
        export_layout.addWidget(self.btn_export_fig)
        
        self.btn_export_pdf = QPushButton("Generer Rapport PDF")
        self.btn_export_pdf.setObjectName("exportBtn")
        self.btn_export_pdf.setEnabled(False)
        self.btn_export_pdf.setCursor(Qt.PointingHandCursor)
        self.btn_export_pdf.clicked.connect(self._export_pdf_report)
        if not REPORTLAB_AVAILABLE:
            self.btn_export_pdf.setToolTip("pip install reportlab")
        export_layout.addWidget(self.btn_export_pdf)
        
        scroll_layout.addWidget(export_group)
        scroll_layout.addStretch()
        
        scroll.setWidget(scroll_widget)
        left_layout.addWidget(scroll)
        
        # =====================================================================
        # PANNEAU DE DROITE
        # =====================================================================
        right_panel = QWidget()
        right_panel.setObjectName("rightPanel")
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(15, 15, 15, 15)
        right_layout.setSpacing(10)
        
        self.tabs = QTabWidget()
        
        # Onglet Visualisation
        viz_tab = QWidget()
        viz_tab_layout = QVBoxLayout(viz_tab)
        viz_tab_layout.setContentsMargins(10, 10, 10, 10)
        
        self.canvas = MatplotlibCanvas(self, width=12, height=9)
        
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setStyleSheet("""
            QToolBar { background: rgba(49, 50, 68, 0.8); border-radius: 8px; padding: 5px; }
            QToolButton { background: transparent; border: none; padding: 5px; }
            QToolButton:hover { background: rgba(137, 180, 250, 0.3); border-radius: 4px; }
        """)
        
        viz_tab_layout.addWidget(self.toolbar)
        viz_tab_layout.addWidget(self.canvas)
        
        self.tabs.addTab(viz_tab, "Visualisation")
        
        # Onglet Statistiques
        stats_tab = QWidget()
        stats_tab_layout = QVBoxLayout(stats_tab)
        stats_tab_layout.setContentsMargins(10, 10, 10, 10)
        
        stats_header = QLabel("Distribution des Metaclusters par Condition")
        stats_header.setStyleSheet("font-size: 12pt; font-weight: 600; color: #89b4fa; padding: 10px;")
        stats_tab_layout.addWidget(stats_header)
        
        self.table_stats = QTableWidget()
        self.table_stats.setColumnCount(4)
        self.table_stats.setHorizontalHeaderLabels(["Metacluster", "NBM (%)", "Patient (%)", "Difference"])
        self.table_stats.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table_stats.setAlternatingRowColors(True)
        stats_tab_layout.addWidget(self.table_stats)
        
        self.tabs.addTab(stats_tab, "Statistiques")
        
        # Onglet Selection Clusters (avec ScrollArea) - Conditionnel selon config
        if self.config.get('onglets', {}).get('clusters_interactifs', True):
            cluster_tab = QWidget()
            cluster_tab_outer_layout = QVBoxLayout(cluster_tab)
            cluster_tab_outer_layout.setContentsMargins(0, 0, 0, 0)
            
            # ScrollArea pour gerer le defilement
            cluster_scroll = QScrollArea()
            cluster_scroll.setWidgetResizable(True)
            cluster_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
            cluster_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
            cluster_scroll.setStyleSheet("""
                QScrollArea {
                    border: none;
                    background: transparent;
                }
            """)
            
            # Widget contenu dans le scroll
            cluster_scroll_widget = QWidget()
            cluster_tab_layout = QVBoxLayout(cluster_scroll_widget)
            cluster_tab_layout.setContentsMargins(10, 10, 10, 10)
            cluster_tab_layout.setSpacing(10)
            
            cluster_header = QLabel("Selection des Clusters a Afficher")
            cluster_header.setStyleSheet("font-size: 12pt; font-weight: 600; color: #89b4fa; padding: 10px;")
            cluster_tab_layout.addWidget(cluster_header)
            
            # Layout horizontal pour liste + star plot
            cluster_main_layout = QHBoxLayout()
            
            # Panneau gauche: liste des clusters
            cluster_left_panel = QWidget()
            cluster_left_layout = QVBoxLayout(cluster_left_panel)
            cluster_left_layout.setContentsMargins(0, 0, 5, 0)
            
            # Boutons selection rapide
            cluster_btn_layout = QHBoxLayout()
            self.btn_select_all = QPushButton("Tout selectionner")
            self.btn_select_all.clicked.connect(self._select_all_clusters)
            cluster_btn_layout.addWidget(self.btn_select_all)
            
            self.btn_deselect_all = QPushButton("Tout deselectionner")
            self.btn_deselect_all.clicked.connect(self._deselect_all_clusters)
            cluster_btn_layout.addWidget(self.btn_deselect_all)
            cluster_left_layout.addLayout(cluster_btn_layout)
            
            # Liste des clusters (selection unique pour le star plot)
            lbl_cluster_select = QLabel("Cliquez sur un cluster pour voir le Star Plot:")
            lbl_cluster_select.setStyleSheet("color: #a6adc8; font-size: 9pt; padding: 5px;")
            cluster_left_layout.addWidget(lbl_cluster_select)
            
            self.cluster_list = QListWidget()
            self.cluster_list.setSelectionMode(QAbstractItemView.SingleSelection)
            self.cluster_list.itemSelectionChanged.connect(self._on_cluster_selection_changed)
            self.cluster_list.setStyleSheet("""
                QListWidget {
                    background: rgba(49, 50, 68, 0.6);
                    border: 1px solid rgba(137, 180, 250, 0.2);
                    border-radius: 8px;
                    padding: 5px;
                }
                QListWidget::item {
                    padding: 8px;
                    border-radius: 4px;
                    margin: 2px;
                }
                QListWidget::item:selected {
                    background: rgba(137, 180, 250, 0.4);
                }
                QListWidget::item:hover {
                    background: rgba(137, 180, 250, 0.2);
                }
            """)
            self.cluster_list.setMaximumWidth(220)
            self.cluster_list.setMinimumHeight(200)
            cluster_left_layout.addWidget(self.cluster_list)
            
            # Info cluster selectionne
            self.lbl_cluster_info = QLabel("Aucun cluster selectionne")
            self.lbl_cluster_info.setStyleSheet("""
                color: #f9e2af;
                font-size: 10pt;
                font-weight: bold;
                background: rgba(249, 226, 175, 0.1);
                padding: 8px;
                border-radius: 6px;
            """)
            self.lbl_cluster_info.setAlignment(Qt.AlignCenter)
            cluster_left_layout.addWidget(self.lbl_cluster_info)
            
            cluster_main_layout.addWidget(cluster_left_panel)
            
            # Panneau droite: Star Plot du cluster
            cluster_right_panel = QWidget()
            cluster_right_layout = QVBoxLayout(cluster_right_panel)
            cluster_right_layout.setContentsMargins(5, 0, 0, 0)
            
            lbl_star_plot = QLabel("Star Plot du Cluster Selectionne")
            lbl_star_plot.setStyleSheet("font-size: 11pt; font-weight: 600; color: #cba6f7; padding: 5px;")
            cluster_right_layout.addWidget(lbl_star_plot)
            
            # Canvas pour le Star Plot
            self.star_canvas = MatplotlibCanvas(self, width=8, height=5)
            self.star_canvas.setMinimumHeight(350)
            cluster_right_layout.addWidget(self.star_canvas)
            
            cluster_main_layout.addWidget(cluster_right_panel, stretch=1)
            
            cluster_tab_layout.addLayout(cluster_main_layout)
            
            # Canvas pour le graphique de comparaison (en bas)
            lbl_comparison = QLabel("Comparaison NBM vs Patient")
            lbl_comparison.setStyleSheet("font-size: 11pt; font-weight: 600; color: #89b4fa; padding: 5px;")
            cluster_tab_layout.addWidget(lbl_comparison)
            
            self.cluster_canvas = MatplotlibCanvas(self, width=10, height=3)
            self.cluster_canvas.setMinimumHeight(250)
            cluster_tab_layout.addWidget(self.cluster_canvas)
            
            # Canvas pour l'image FlowSOM Stars MST (en bas)
            fsom_stars_header = QHBoxLayout()
            lbl_fsom_stars = QLabel("FlowSOM MST - Star Plot du Cluster")
            lbl_fsom_stars.setStyleSheet("font-size: 11pt; font-weight: 600; color: #a6e3a1; padding: 5px;")
            fsom_stars_header.addWidget(lbl_fsom_stars)
            
            self.btn_generate_mst = QPushButton("Generer MST")
            self.btn_generate_mst.setObjectName("successBtn")
            self.btn_generate_mst.setFixedWidth(150)
            self.btn_generate_mst.setEnabled(False)
            self.btn_generate_mst.clicked.connect(self._update_fsom_stars_mst)
            fsom_stars_header.addWidget(self.btn_generate_mst)
            fsom_stars_header.addStretch()
            
            cluster_tab_layout.addLayout(fsom_stars_header)
            
            self.fsom_stars_canvas = MatplotlibCanvas(self, width=10, height=6)
            self.fsom_stars_canvas.setMinimumHeight(450)
            cluster_tab_layout.addWidget(self.fsom_stars_canvas)
            
            # Ajouter un stretch a la fin pour eviter l'etirement
            cluster_tab_layout.addStretch()
            
            # Finaliser le scroll
            cluster_scroll.setWidget(cluster_scroll_widget)
            cluster_tab_outer_layout.addWidget(cluster_scroll)
            
            self.tabs.addTab(cluster_tab, "Clusters Interactifs")
        else:
            # Widgets par defaut si onglet desactive
            self.cluster_list = None
            self.lbl_cluster_info = None
            self.star_canvas = None
            self.cluster_canvas = None
            self.fsom_stars_canvas = None
            self.btn_select_all = None
            self.btn_deselect_all = None
            self.btn_generate_mst = None
        
        # Onglet Visualisation FCS - Style Kaluza
        if self.config.get('onglets', {}).get('fcs_viewer', True):
            fcs_tab = QWidget()
            fcs_tab_layout = QVBoxLayout(fcs_tab)
            fcs_tab_layout.setContentsMargins(10, 10, 10, 10)
            
            fcs_header = QLabel("Visualisation FCS - Cytometrie")
            fcs_header.setStyleSheet("font-size: 12pt; font-weight: 600; color: #89b4fa; padding: 10px;")
            fcs_tab_layout.addWidget(fcs_header)
            
            # Panel de controle
            fcs_control_panel = QWidget()
            fcs_control_layout = QHBoxLayout(fcs_control_panel)
            fcs_control_layout.setContentsMargins(0, 0, 0, 0)
            
            # Bouton charger FCS
            self.btn_load_fcs_viz = QPushButton("Charger FCS")
            self.btn_load_fcs_viz.setObjectName("primaryBtn")
            self.btn_load_fcs_viz.setCursor(Qt.PointingHandCursor)
            self.btn_load_fcs_viz.clicked.connect(self._load_fcs_for_visualization)
            fcs_control_layout.addWidget(self.btn_load_fcs_viz)
            
            # Selection axe X
            lbl_x = QLabel("Axe X:")
            lbl_x.setStyleSheet("color: #cdd6f4; font-weight: 600;")
            fcs_control_layout.addWidget(lbl_x)
            
            self.combo_fcs_x = QComboBox()
            self.combo_fcs_x.setMinimumWidth(120)
            self.combo_fcs_x.currentIndexChanged.connect(self._update_fcs_plot)
            fcs_control_layout.addWidget(self.combo_fcs_x)
            
            # Selection axe Y
            lbl_y = QLabel("Axe Y:")
            lbl_y.setStyleSheet("color: #cdd6f4; font-weight: 600;")
            fcs_control_layout.addWidget(lbl_y)
            
            self.combo_fcs_y = QComboBox()
            self.combo_fcs_y.setMinimumWidth(120)
            self.combo_fcs_y.currentIndexChanged.connect(self._update_fcs_plot)
            fcs_control_layout.addWidget(self.combo_fcs_y)
            
            # Type de plot
            lbl_plot = QLabel("Type:")
            lbl_plot.setStyleSheet("color: #cdd6f4; font-weight: 600;")
            fcs_control_layout.addWidget(lbl_plot)
            
            self.combo_fcs_plot_type = QComboBox()
            self.combo_fcs_plot_type.addItems(["Scatter", "Densite", "Contour"])
            self.combo_fcs_plot_type.currentIndexChanged.connect(self._update_fcs_plot)
            fcs_control_layout.addWidget(self.combo_fcs_plot_type)
            
            # Coloration par cluster/metacluster
            lbl_color = QLabel("Couleur:")
            lbl_color.setStyleSheet("color: #cdd6f4; font-weight: 600;")
            fcs_control_layout.addWidget(lbl_color)
            
            self.combo_fcs_color = QComboBox()
            self.combo_fcs_color.addItems(["Aucune", "FlowSOM_cluster", "FlowSOM_metacluster", "Condition"])
            self.combo_fcs_color.setToolTip("Colorier les points par cluster ou metacluster")
            self.combo_fcs_color.currentIndexChanged.connect(self._update_fcs_plot)
            fcs_control_layout.addWidget(self.combo_fcs_color)
            
            # Sous-echantillonnage
            lbl_sample = QLabel("Cellules:")
            lbl_sample.setStyleSheet("color: #cdd6f4; font-weight: 600;")
            fcs_control_layout.addWidget(lbl_sample)
            
            self.spin_fcs_cells = QSpinBox()
            self.spin_fcs_cells.setRange(1000, 500000)
            self.spin_fcs_cells.setValue(10000)
            self.spin_fcs_cells.setSingleStep(5000)
            fcs_control_layout.addWidget(self.spin_fcs_cells)
            
            # Checkbox toutes les cellules
            self.chk_fcs_all_cells = QCheckBox("Toutes")
            self.chk_fcs_all_cells.setToolTip("Afficher toutes les cellules (peut etre lent)")
            self.chk_fcs_all_cells.stateChanged.connect(self._toggle_fcs_all_cells)
            fcs_control_layout.addWidget(self.chk_fcs_all_cells)
            
            # Checkbox jitter (dispersion) pour coordonnees SOM
            self.chk_fcs_jitter = QCheckBox("Jitter")
            self.chk_fcs_jitter.setToolTip("Ajouter une dispersion aux coordonnees SOM (xGrid, yGrid, xNodes, yNodes)")
            self.chk_fcs_jitter.setChecked(True)  # Active par defaut
            self.chk_fcs_jitter.stateChanged.connect(self._update_fcs_plot)
            fcs_control_layout.addWidget(self.chk_fcs_jitter)
            
            # Bouton rafraichir
            self.btn_refresh_fcs = QPushButton("Rafraichir")
            self.btn_refresh_fcs.setObjectName("secondaryBtn")
            self.btn_refresh_fcs.setCursor(Qt.PointingHandCursor)
            self.btn_refresh_fcs.clicked.connect(self._update_fcs_plot)
            fcs_control_layout.addWidget(self.btn_refresh_fcs)
            
            fcs_control_layout.addStretch()
            fcs_tab_layout.addWidget(fcs_control_panel)
            
            # Canvas pour le plot FCS
            self.fcs_viz_canvas = MatplotlibCanvas(self, width=10, height=8, dpi=100)
            self.fcs_viz_canvas.setMinimumHeight(500)
            fcs_tab_layout.addWidget(self.fcs_viz_canvas)
            
            # Label d'info
            self.lbl_fcs_info = QLabel("Chargez un fichier FCS pour visualiser")
            self.lbl_fcs_info.setStyleSheet("color: #a6adc8; padding: 5px;")
            fcs_tab_layout.addWidget(self.lbl_fcs_info)
            
            self.tabs.addTab(fcs_tab, "Visualisation FCS")
        else:
            self.fcs_viz_canvas = None
            self.combo_fcs_x = None
            self.combo_fcs_y = None
        
        # Onglet Suivi Patient - Conditionnel selon config
        if self.config.get('onglets', {}).get('suivi_patient', True):
            patient_tab = QWidget()
            patient_tab_layout = QVBoxLayout(patient_tab)
            patient_tab_layout.setContentsMargins(10, 10, 10, 10)
            
            patient_header = QLabel("Suivi Longitudinal Patient")
            patient_header.setStyleSheet("font-size: 12pt; font-weight: 600; color: #89b4fa; padding: 10px;")
            patient_tab_layout.addWidget(patient_header)
            
            # Panel de controle patient
            patient_control = QWidget()
            patient_control_layout = QHBoxLayout(patient_control)
            patient_control_layout.setContentsMargins(0, 0, 0, 0)
            
            # ID Patient
            lbl_patient_id = QLabel("ID Patient:")
            lbl_patient_id.setStyleSheet("color: #cdd6f4; font-weight: 600;")
            patient_control_layout.addWidget(lbl_patient_id)
            
            self.edit_patient_id = QLineEdit()
            self.edit_patient_id.setPlaceholderText("Ex: PAT001")
            self.edit_patient_id.setMaximumWidth(120)
            patient_control_layout.addWidget(self.edit_patient_id)
            
            # Boutons sauvegarder/charger
            self.btn_save_patient = QPushButton("Sauvegarder Patient")
            self.btn_save_patient.setObjectName("successBtn")
            self.btn_save_patient.setCursor(Qt.PointingHandCursor)
            self.btn_save_patient.clicked.connect(self._save_patient_data)
            patient_control_layout.addWidget(self.btn_save_patient)
            
            self.btn_load_patient = QPushButton("Charger Patient")
            self.btn_load_patient.setObjectName("primaryBtn")
            self.btn_load_patient.setCursor(Qt.PointingHandCursor)
            self.btn_load_patient.clicked.connect(self._load_patient_data)
            patient_control_layout.addWidget(self.btn_load_patient)
            
            # Ajouter timepoint
            self.btn_add_timepoint = QPushButton("Ajouter Timepoint")
            self.btn_add_timepoint.setObjectName("secondaryBtn")
            self.btn_add_timepoint.setCursor(Qt.PointingHandCursor)
            self.btn_add_timepoint.clicked.connect(self._add_patient_timepoint)
            self.btn_add_timepoint.setEnabled(False)
            patient_control_layout.addWidget(self.btn_add_timepoint)
            
            patient_control_layout.addStretch()
            patient_tab_layout.addWidget(patient_control)
            
            # Liste des timepoints
            self.patient_timeline = QTableWidget()
            self.patient_timeline.setColumnCount(5)
            self.patient_timeline.setHorizontalHeaderLabels(["Timepoint", "Date", "Fichier", "Clusters", "Notes"])
            self.patient_timeline.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            self.patient_timeline.setAlternatingRowColors(True)
            patient_tab_layout.addWidget(self.patient_timeline)
            
            # Canvas pour l'evolution temporelle
            self.patient_canvas = MatplotlibCanvas(self, width=10, height=5)
            patient_tab_layout.addWidget(self.patient_canvas)
            
            self.tabs.addTab(patient_tab, "Suivi Patient")
        else:
            self.patient_timeline = None
            self.patient_canvas = None
            self.edit_patient_id = None
        
        # Onglet Logs - Conditionnel selon config
        if self.config.get('onglets', {}).get('logs', True):
            logs_tab = QWidget()
            logs_tab_layout = QVBoxLayout(logs_tab)
            logs_tab_layout.setContentsMargins(10, 10, 10, 10)
            
            logs_header = QLabel("Journal d'Analyse")
            logs_header.setStyleSheet("font-size: 12pt; font-weight: 600; color: #89b4fa; padding: 10px;")
            logs_tab_layout.addWidget(logs_header)
            
            self.text_logs = QTextEdit()
            self.text_logs.setReadOnly(True)
            self.text_logs.setPlaceholderText("Les logs apparaitront ici...")
            logs_tab_layout.addWidget(self.text_logs)
            
            self.tabs.addTab(logs_tab, "Logs")
        else:
            self.text_logs = None
        
        right_layout.addWidget(self.tabs)
        
        # Ajout au splitter
        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        
    def _toggle_auto_cluster(self, state):
        self.spin_n_clusters.setEnabled(state != Qt.Checked)
    
    def _toggle_nbm_mode(self, state):
        """Active/desactive le mode de comparaison avec NBM."""
        compare_mode = state == Qt.Checked
        self.btn_healthy.setEnabled(compare_mode)
        
        if not compare_mode:
            # Effacer la reference NBM si on passe en mode patient seul
            self.healthy_files = []
            self.healthy_folder = None
            self.lbl_healthy.setText("Mode: Patient seul (sans reference)")
            self.lbl_healthy.setStyleSheet("color: #a6adc8; font-style: italic;")
        else:
            self.lbl_healthy.setText("Aucune reference selectionnee")
            self.lbl_healthy.setStyleSheet("")
        
        self._update_file_info()
        
    def _toggle_all_cells(self, state):
        """Active/desactive le spinbox du nombre de cellules selon l'option 'toutes les cellules'."""
        self.spin_embedding_cells.setEnabled(state != Qt.Checked)
        if state == Qt.Checked:
            self._log("UMAP/t-SNE: mode toutes les cellules active")
    
    def _toggle_fcs_all_cells(self, state):
        """Active/desactive le spinbox du nombre de cellules FCS."""
        if hasattr(self, 'spin_fcs_cells'):
            self.spin_fcs_cells.setEnabled(state != Qt.Checked)
        # Rafraichir le plot automatiquement
        self._update_fcs_plot()
        
    def _load_healthy_folder(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Selectionner Reference NBM (Moelle Normale)", "",
            QFileDialog.ShowDirsOnly
        )
        if folder:
            self.healthy_folder = folder
            self.healthy_files = self._get_fcs_files(folder)
            n = len(self.healthy_files)
            self.lbl_healthy.setText(f"{Path(folder).name} ({n} fichiers)")
            self.lbl_healthy.setStyleSheet("color: #a6e3a1; font-weight: 600; padding: 5px 10px; background: rgba(166, 227, 161, 0.1); border-radius: 6px;")
            self._update_file_info()
            self._log(f"Reference NBM: {folder} ({n} fichiers)")
            
    def _load_pathological_folder(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Selectionner Patient Diagnostic", "",
            QFileDialog.ShowDirsOnly
        )
        if folder:
            self.pathological_folder = folder
            self.pathological_files = self._get_fcs_files(folder)
            n = len(self.pathological_files)
            self.lbl_pathological.setText(f"{Path(folder).name} ({n} fichiers)")
            self.lbl_pathological.setStyleSheet("color: #f38ba8; font-weight: 600; padding: 5px 10px; background: rgba(243, 139, 168, 0.1); border-radius: 6px;")
            self._update_file_info()
            self._log(f"Patient diagnostic: {folder} ({n} fichiers)")
            
    def _get_fcs_files(self, folder):
        folder_path = Path(folder)
        # Utiliser un set pour éviter les doublons (Windows est insensible à la casse)
        files = set()
        for f in folder_path.glob("*.fcs"):
            files.add(str(f))
        for f in folder_path.glob("*.FCS"):
            files.add(str(f))
        return list(files)
        
    def _update_file_info(self):
        n_h = len(self.healthy_files)
        n_p = len(self.pathological_files)
        
        # Verifier le mode d'analyse
        compare_mode = hasattr(self, 'chk_compare_nbm') and self.chk_compare_nbm.isChecked()
        
        if compare_mode:
            self.lbl_file_info.setText(f"Fichiers: {n_h} reference(s), {n_p} patient(s)")
            # En mode comparaison, il faut les deux
            self.btn_run.setEnabled(n_h > 0 and n_p > 0)
        else:
            self.lbl_file_info.setText(f"Fichiers: {n_p} patient(s) (mode patient seul)")
            # En mode patient seul, juste le patient suffit
            self.btn_run.setEnabled(n_p > 0)
    
    # =========================================================================
    # VISUALISATION FCS (style Kaluza)
    # =========================================================================
    
    def _read_fcs_binary(self, file_path: str):
        """
        Lecture binaire directe d'un fichier FCS.
        Fallback robuste quand les autres méthodes échouent.
        """
        import struct
        
        with open(file_path, 'rb') as f:
            # Lire le header FCS (58 bytes)
            header = f.read(58)
            
            # Version FCS
            version = header[0:6].decode('ascii').strip()
            
            # Positions du segment TEXT
            text_start = int(header[10:18].decode('ascii').strip())
            text_end = int(header[18:26].decode('ascii').strip())
            
            # Positions du segment DATA
            data_start = int(header[26:34].decode('ascii').strip())
            data_end = int(header[34:42].decode('ascii').strip())
            
            # Lire le segment TEXT
            f.seek(text_start)
            text_segment = f.read(text_end - text_start + 1)
            
            # Decoder le segment TEXT (latin-1 pour FCS)
            try:
                text_str = text_segment.decode('latin-1')
            except:
                text_str = text_segment.decode('utf-8', errors='replace')
            
            # Parser les paires cle-valeur
            delimiter = text_str[0]
            parts = text_str[1:].split(delimiter)
            text_dict = {}
            for i in range(0, len(parts) - 1, 2):
                key = parts[i].strip().upper()
                value = parts[i + 1].strip() if i + 1 < len(parts) else ""
                text_dict[key] = value
            
            # Obtenir les parametres essentiels
            n_params = int(text_dict.get('$PAR', text_dict.get('PAR', 0)))
            n_events = int(text_dict.get('$TOT', text_dict.get('TOT', 0)))
            datatype = text_dict.get('$DATATYPE', text_dict.get('DATATYPE', 'F')).upper()
            byteord = text_dict.get('$BYTEORD', text_dict.get('BYTEORD', '1,2,3,4'))
            
            if n_params == 0 or n_events == 0:
                raise ValueError(f"Parametres invalides: {n_params} params, {n_events} events")
            
            # Determiner l'endianness
            if byteord in ['1,2,3,4', '1,2']:
                endian = '<'  # Little endian
            else:
                endian = '>'  # Big endian
            
            # Noms des canaux
            channel_names = []
            for i in range(1, n_params + 1):
                # Chercher le nom ($PnS ou $PnN)
                name = None
                for key in [f'$P{i}S', f'P{i}S', f'$P{i}N', f'P{i}N']:
                    if key in text_dict:
                        name = text_dict[key]
                        break
                if not name:
                    name = f'Channel_{i}'
                channel_names.append(name)
            
            # Lire les donnees
            f.seek(data_start)
            data_bytes = f.read(data_end - data_start + 1)
            
            # Decoder selon le type de donnees
            if datatype == 'F':
                # Float 32 bits
                fmt = f'{endian}{n_params}f'
                bytes_per_event = n_params * 4
            elif datatype == 'D':
                # Double 64 bits
                fmt = f'{endian}{n_params}d'
                bytes_per_event = n_params * 8
            elif datatype == 'I':
                # Integer - determiner la taille par $PnB
                bits = int(text_dict.get('$P1B', text_dict.get('P1B', 16)))
                if bits == 16:
                    fmt = f'{endian}{n_params}H'
                    bytes_per_event = n_params * 2
                else:  # 32 bits
                    fmt = f'{endian}{n_params}I'
                    bytes_per_event = n_params * 4
            else:
                raise ValueError(f"Type de donnees non supporte: {datatype}")
            
            # Lire tous les events
            events = []
            for i in range(n_events):
                offset = i * bytes_per_event
                if offset + bytes_per_event <= len(data_bytes):
                    try:
                        event = struct.unpack(fmt, data_bytes[offset:offset + bytes_per_event])
                        events.append(event)
                    except:
                        break
            
            if len(events) == 0:
                raise ValueError("Aucun event lu")
            
            # Creer l'AnnData
            data_array = np.array(events, dtype=np.float32)
            adata = ad.AnnData(data_array)
            adata.var_names = channel_names
            
            self._log(f"Lecture binaire: {len(events)} events, {n_params} canaux")
            return adata
    
    def _load_fcs_for_visualization(self):
        """Charge un fichier FCS pour la visualisation interactive."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Charger un fichier FCS", "",
            "FCS Files (*.fcs *.FCS)"
        )
        if not file_path:
            return
            
        try:
            self._log(f"Chargement FCS pour visualisation: {Path(file_path).name}")
            
            # Essayer plusieurs methodes de chargement
            adata = None
            last_error = None
            
            # Methode 1: flowsom (fichiers FCS standards)
            try:
                adata = fs.io.read_FCS(file_path)
                self._log("Charge avec flowsom")
            except Exception as e1:
                self._log(f"Methode flowsom echouee: {str(e1)[:50]}")
                last_error = e1
            
            # Methode 2: flowio (tres permissif, avec options relaxees)
            if adata is None:
                try:
                    import flowio
                    # Les fichiers FCS utilisent l'encodage latin-1, pas UTF-8
                    fcs_data = flowio.FlowData(file_path)
                    
                    # Recuperer les events
                    events = np.reshape(fcs_data.events, (-1, fcs_data.channel_count))
                    n_events, n_channels = events.shape
                    
                    if n_events > 0 and n_channels > 0:
                        # Recuperer les noms des canaux avec toutes les variantes possibles
                        channel_names = []
                        for i in range(1, n_channels + 1):
                            # Essayer toutes les variantes de cles (avec et sans $)
                            possible_keys = [
                                f'$P{i}N', f'P{i}N', 
                                f'$P{i}S', f'P{i}S',
                                f'p{i}n', f'p{i}s'
                            ]
                            name = None
                            for key in possible_keys:
                                if key in fcs_data.text:
                                    name = fcs_data.text[key]
                                    break
                            if name is None:
                                name = f'Channel_{i}'
                            channel_names.append(str(name))
                        
                        adata = ad.AnnData(events.astype(np.float32))
                        adata.var_names = channel_names
                        self._log("Charge avec flowio")
                    else:
                        raise ValueError(f"Fichier vide: {n_events} events, {n_channels} channels")
                except ImportError:
                    self._log("flowio non installe")
                except Exception as e2:
                    self._log(f"Methode flowio echouee: {str(e2)[:50]}")
                    last_error = e2
            
            # Methode 3: fcsparser (plus permissif)
            if adata is None:
                try:
                    import fcsparser
                    meta, data = fcsparser.parse(file_path, meta_data_only=False, 
                                                  reformat_meta=False, 
                                                  channel_naming='$PnS')
                    adata = ad.AnnData(data.values.astype(np.float32))
                    adata.var_names = list(data.columns)
                    self._log("Charge avec fcsparser")
                except ImportError:
                    self._log("fcsparser non installe")
                except Exception as e3:
                    # Essayer avec un autre channel_naming
                    try:
                        import fcsparser
                        meta, data = fcsparser.parse(file_path, meta_data_only=False, 
                                                      reformat_meta=False, 
                                                      channel_naming='$PnN')
                        adata = ad.AnnData(data.values.astype(np.float32))
                        adata.var_names = list(data.columns)
                        self._log("Charge avec fcsparser (PnN)")
                    except:
                        self._log(f"Methode fcsparser echouee: {str(e3)[:50]}")
                        last_error = e3
            
            # Methode 4: Lecture binaire directe (fallback ultime)
            if adata is None:
                try:
                    adata = self._read_fcs_binary(file_path)
                    if adata is not None:
                        self._log("Charge avec lecture binaire directe")
                except Exception as e4:
                    self._log(f"Methode binaire echouee: {str(e4)[:50]}")
                    last_error = e4
            
            if adata is None:
                raise Exception(f"Impossible de charger le FCS. Derniere erreur: {last_error}")
            
            self.current_fcs_adata = adata
            
            # Remplir les combos avec les marqueurs
            markers = list(adata.var_names)
            
            self.combo_fcs_x.blockSignals(True)
            self.combo_fcs_y.blockSignals(True)
            
            self.combo_fcs_x.clear()
            self.combo_fcs_y.clear()
            self.combo_fcs_x.addItems(markers)
            self.combo_fcs_y.addItems(markers)
            
            # Mettre a jour le combo de coloration avec les colonnes de clustering disponibles
            if hasattr(self, 'combo_fcs_color'):
                self.combo_fcs_color.blockSignals(True)
                self.combo_fcs_color.clear()
                self.combo_fcs_color.addItem("Aucune")
                
                # Ajouter les colonnes de clustering/metaclustering si presentes
                # Recherche insensible a la casse
                color_patterns = ['flowsom_cluster', 'flowsom_metacluster', 'condition', 
                                  'cluster', 'metacluster', 'flowsom']
                markers_lower = [m.lower() for m in markers]
                
                for marker in markers:
                    marker_lower = marker.lower()
                    # Verifier si c'est une colonne de clustering/coloration
                    if any(pattern in marker_lower for pattern in color_patterns):
                        self.combo_fcs_color.addItem(marker)
                
                # Log pour debug
                n_color_options = self.combo_fcs_color.count()
                self._log(f"Options de coloration disponibles: {n_color_options - 1}")
                for i in range(1, n_color_options):
                    self._log(f"   - {self.combo_fcs_color.itemText(i)}")
                
                self.combo_fcs_color.blockSignals(False)
            
            # Selectionner FSC-A et SSC-A par defaut si disponibles
            fsc_idx = next((i for i, m in enumerate(markers) if 'FSC' in m.upper()), 0)
            ssc_idx = next((i for i, m in enumerate(markers) if 'SSC' in m.upper()), min(1, len(markers)-1))
            
            self.combo_fcs_x.setCurrentIndex(fsc_idx)
            self.combo_fcs_y.setCurrentIndex(ssc_idx)
            
            self.combo_fcs_x.blockSignals(False)
            self.combo_fcs_y.blockSignals(False)
            
            n_cells = adata.shape[0]
            n_markers = adata.shape[1]
            self.lbl_fcs_info.setText(f"Fichier: {Path(file_path).name} | {n_cells:,} cellules | {n_markers} parametres")
            
            self._update_fcs_plot()
            self._log(f"FCS charge: {n_cells:,} cellules, {n_markers} parametres")
            
        except Exception as e:
            QMessageBox.critical(self, "Erreur", f"Erreur de chargement FCS:\n{str(e)}")
            self._log(f"Erreur chargement FCS: {str(e)}")
    
    def _update_fcs_plot(self):
        """Met a jour le plot FCS avec les parametres selectionnes."""
        if self.current_fcs_adata is None:
            return
            
        if self.fcs_viz_canvas is None:
            return
        
        try:
            x_marker = self.combo_fcs_x.currentText()
            y_marker = self.combo_fcs_y.currentText()
            plot_type = self.combo_fcs_plot_type.currentText()
            
            # Coloration par cluster/metacluster
            color_by = "Aucune"
            if hasattr(self, 'combo_fcs_color'):
                color_by = self.combo_fcs_color.currentText()
            
            # Determiner si on affiche toutes les cellules
            show_all = hasattr(self, 'chk_fcs_all_cells') and self.chk_fcs_all_cells.isChecked()
            max_cells = self.spin_fcs_cells.value() if not show_all else float('inf')
            
            if not x_marker or not y_marker:
                return
            
            # Extraire les donnees
            X = self.current_fcs_adata.X
            if hasattr(X, 'toarray'):
                X = X.toarray()
            
            var_names = list(self.current_fcs_adata.var_names)
            x_idx = var_names.index(x_marker)
            y_idx = var_names.index(y_marker)
            
            x_data = X[:, x_idx]
            y_data = X[:, y_idx]
            
            # Preparer les donnees de coloration si demande
            color_data = None
            if color_by != "Aucune" and color_by in var_names:
                color_idx = var_names.index(color_by)
                color_data = X[:, color_idx]
            
            # Detecter si ce sont des coordonnees SOM (qui necessitent du jitter)
            som_grid_cols = ['xgrid', 'ygrid']  # Coordonnees grille (entieres)
            som_nodes_cols = ['xnodes', 'ynodes']  # Coordonnees MST (continues)
            is_grid_x = x_marker.lower() in som_grid_cols
            is_grid_y = y_marker.lower() in som_grid_cols
            is_nodes_x = x_marker.lower() in som_nodes_cols
            is_nodes_y = y_marker.lower() in som_nodes_cols
            is_som_x = is_grid_x or is_nodes_x
            is_som_y = is_grid_y or is_nodes_y
            
            # Appliquer le jitter si demande et si ce sont des coordonnees SOM
            apply_jitter = hasattr(self, 'chk_fcs_jitter') and self.chk_fcs_jitter.isChecked()
            if apply_jitter and (is_som_x or is_som_y):
                # Jitter circulaire (polaire) pour former des disques comme Kaluza
                n_points = len(x_data)
                # Rayon aleatoire avec distribution sqrt pour repartition uniforme dans le cercle
                r = np.sqrt(np.random.uniform(0, 1, n_points))
                # Angle aleatoire
                theta = np.random.uniform(0, 2 * np.pi, n_points)
                
                # Amplitude du jitter selon le type de coordonnees (style Kaluza)
                if is_grid_x or is_grid_y:
                    # Pour la grille (valeurs entieres 0-10): rayon ~0.35 pour cercles bien separes
                    jitter_radius = 0.35
                else:
                    # Pour les coordonnees MST (xNodes/yNodes echelle ~1000): rayon fixe petit
                    # Kaluza utilise un jitter d'environ 20-25 unites sur echelle 1000
                    jitter_radius = 20.0
                
                # Appliquer le jitter circulaire
                if is_som_x:
                    x_data = x_data + r * np.cos(theta) * jitter_radius
                if is_som_y:
                    y_data = y_data + r * np.sin(theta) * jitter_radius
            
            # Filtrer les valeurs manquantes (-999) pour les colonnes de reduction dimensionnelle
            # Cela evite d'afficher les cellules sans coordonnees valides
            dim_reduction_cols = ['tSNE1', 'tSNE2', 'UMAP1', 'UMAP2', 'tsne1', 'tsne2', 'umap1', 'umap2']
            is_dim_reduction_x = x_marker.lower() in [c.lower() for c in dim_reduction_cols]
            is_dim_reduction_y = y_marker.lower() in [c.lower() for c in dim_reduction_cols]
            
            # Masque de validite: exclure -999, NaN, Inf
            MISSING_VALUE = -999.0
            valid_x = np.isfinite(x_data) & (x_data != MISSING_VALUE)
            valid_y = np.isfinite(y_data) & (y_data != MISSING_VALUE)
            
            if is_dim_reduction_x or is_dim_reduction_y:
                # Pour les colonnes de reduction, filtrer strictement les valeurs manquantes
                initial_mask = valid_x & valid_y
            else:
                # Pour les autres colonnes, juste filtrer NaN/Inf
                initial_mask = np.isfinite(x_data) & np.isfinite(y_data)
            
            # Appliquer le masque a toutes les donnees
            x_data = x_data[initial_mask]
            y_data = y_data[initial_mask]
            if color_data is not None:
                color_data = color_data[initial_mask]
            n_cells_valid = len(x_data)
            
            if is_dim_reduction_x or is_dim_reduction_y:
                self._log(f"FCS plot: {n_cells_valid:,} cellules avec coordonnees valides sur {len(X):,}")
            
            # Sous-echantillonnage (sauf si toutes les cellules)
            n_cells = len(x_data)
            if not show_all and n_cells > max_cells:
                idx = np.random.choice(n_cells, int(max_cells), replace=False)
                x_data = x_data[idx]
                y_data = y_data[idx]
                if color_data is not None:
                    color_data = color_data[idx]
            
            # Mettre a jour l'info
            n_total_in_file = len(X)
            if hasattr(self, 'lbl_fcs_info') and self.lbl_fcs_info is not None:
                if is_dim_reduction_x or is_dim_reduction_y:
                    self.lbl_fcs_info.setText(f"Affichage: {len(x_data):,} / {n_cells_valid:,} valides (total: {n_total_in_file:,})")
                else:
                    self.lbl_fcs_info.setText(f"Affichage: {len(x_data):,} / {n_total_in_file:,} cellules")
            
            # Effacer le canvas
            self.fcs_viz_canvas.clear_figure()
            ax = self.fcs_viz_canvas.fig.add_subplot(111)
            
            # Preparer la palette de couleurs pour clusters/metaclusters
            scatter_colors = '#89b4fa'  # Couleur par defaut
            legend_handles = None
            
            if color_data is not None and plot_type == "Scatter":
                # Coloration par cluster/metacluster
                unique_values = np.unique(color_data[np.isfinite(color_data)])
                n_colors = len(unique_values)
                
                # Choisir la colormap selon le nombre de valeurs
                if n_colors <= 20:
                    cmap = plt.cm.tab20
                elif n_colors <= 40:
                    cmap = plt.cm.tab20b
                else:
                    cmap = plt.cm.turbo
                
                # Mapper les valeurs a des couleurs
                color_indices = np.searchsorted(unique_values, color_data)
                scatter_colors = cmap(color_indices / max(n_colors - 1, 1))
                
                # Creer la legende (limiter a 20 elements)
                from matplotlib.patches import Patch
                if n_colors <= 20:
                    legend_handles = []
                    for i, val in enumerate(unique_values):
                        label = f"{color_by.replace('FlowSOM_', '')} {int(val)}"
                        legend_handles.append(Patch(facecolor=cmap(i / max(n_colors - 1, 1)), 
                                                    edgecolor='white', label=label))
            
            # Tracer selon le type
            if plot_type == "Scatter":
                scatter = ax.scatter(x_data, y_data, s=3, alpha=0.6, c=scatter_colors, 
                          edgecolors='none', rasterized=True)
                
                # Ajouter la legende si coloration par cluster
                if legend_handles is not None:
                    ax.legend(handles=legend_handles, loc='upper right', 
                             fontsize=7, facecolor='#313244', labelcolor='#cdd6f4',
                             edgecolor='#45475a', framealpha=0.9,
                             ncol=2 if len(legend_handles) > 10 else 1)
                             
            elif plot_type == "Densite":
                # Histogramme 2D (density plot)
                h = ax.hist2d(x_data, y_data, bins=100, cmap='viridis', 
                             norm=plt.matplotlib.colors.LogNorm())
                cbar = self.fcs_viz_canvas.fig.colorbar(h[3], ax=ax, label='Densité')
                cbar.ax.tick_params(colors='#cdd6f4')
                cbar.ax.yaxis.label.set_color('#cdd6f4')
            elif plot_type == "Contour":
                # Contour plot
                from scipy import stats
                try:
                    # Estimation densite par noyau
                    xmin, xmax = x_data.min(), x_data.max()
                    ymin, ymax = y_data.min(), y_data.max()
                    
                    # Sous-echantillonner pour le calcul KDE (plus rapide)
                    n_kde = min(5000, len(x_data))
                    kde_idx = np.random.choice(len(x_data), n_kde, replace=False)
                    
                    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
                    positions = np.vstack([xx.ravel(), yy.ravel()])
                    values = np.vstack([x_data[kde_idx], y_data[kde_idx]])
                    kernel = stats.gaussian_kde(values)
                    f = np.reshape(kernel(positions).T, xx.shape)
                    
                    ax.contourf(xx, yy, f, levels=20, cmap='viridis')
                    ax.contour(xx, yy, f, levels=10, colors='white', linewidths=0.3, alpha=0.5)
                except Exception:
                    # Fallback sur scatter si erreur KDE
                    ax.scatter(x_data, y_data, s=2, alpha=0.5, c='#89b4fa', 
                              edgecolors='none', rasterized=True)
            
            # Ajuster les limites avec marges
            x_margin = (x_data.max() - x_data.min()) * 0.02
            y_margin = (y_data.max() - y_data.min()) * 0.02
            ax.set_xlim(x_data.min() - x_margin, x_data.max() + x_margin)
            ax.set_ylim(y_data.min() - y_margin, y_data.max() + y_margin)
            
            # Labels et titre
            ax.set_xlabel(x_marker, color='#cdd6f4', fontsize=12, fontweight='bold')
            ax.set_ylabel(y_marker, color='#cdd6f4', fontsize=12, fontweight='bold')
            
            # Construire le titre avec info sur jitter et coloration
            title_parts = [f"{x_marker} vs {y_marker}"]
            subtitle_parts = [f"{len(x_data):,} cellules"]
            if apply_jitter and (is_som_x or is_som_y):
                subtitle_parts.append("jitter")
            if color_by != "Aucune":
                subtitle_parts.append(f"couleur: {color_by.replace('FlowSOM_', '')}")
            
            ax.set_title(f"{title_parts[0]}\n{' | '.join(subtitle_parts)}", 
                        fontsize=13, color='#cdd6f4', fontweight='bold', pad=15)
            
            # Formatage des axes avec notation scientifique si nécessaire
            from matplotlib.ticker import ScalarFormatter, FuncFormatter
            
            def format_axis(value, pos):
                if abs(value) >= 1e6:
                    return f'{value/1e6:.1f}M'
                elif abs(value) >= 1e3:
                    return f'{value/1e3:.0f}K'
                else:
                    return f'{value:.0f}'
            
            ax.xaxis.set_major_formatter(FuncFormatter(format_axis))
            ax.yaxis.set_major_formatter(FuncFormatter(format_axis))
            
            # Style des axes
            ax.tick_params(colors='#cdd6f4', labelsize=10)
            ax.set_facecolor('#1e1e2e')
            ax.spines['bottom'].set_color('#45475a')
            ax.spines['left'].set_color('#45475a')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            # Grille légère
            ax.grid(True, alpha=0.15, color='#6c7086', linestyle='--', linewidth=0.5)
            
            self.fcs_viz_canvas.fig.patch.set_facecolor('#1e1e2e')
            self.fcs_viz_canvas.fig.tight_layout(pad=1.5)
            self.fcs_viz_canvas.draw()
            
        except Exception as e:
            self._log(f"Erreur plot FCS: {str(e)}")
    
    # =========================================================================
    # SUIVI PATIENT - SAUVEGARDE / CHARGEMENT
    # =========================================================================
    
    def _get_patient_data_folder(self) -> Path:
        """Retourne le dossier de stockage des donnees patient."""
        # Creer le dossier Patients dans le dossier de l'application
        data_folder = Path(__file__).parent / "Patients"
        data_folder.mkdir(exist_ok=True)
        return data_folder
    
    def _save_patient_data(self):
        """Sauvegarde les donnees du patient actuel."""
        if self.edit_patient_id is None:
            return
            
        patient_id = self.edit_patient_id.text().strip()
        if not patient_id:
            QMessageBox.warning(self, "Attention", "Entrez un ID patient avant de sauvegarder.")
            return
        
        if not self.result:
            QMessageBox.warning(self, "Attention", "Lancez d'abord une analyse FlowSOM.")
            return
        
        try:
            # Preparer les donnees patient
            patient_data = {
                'patient_id': patient_id,
                'created_at': datetime.now().isoformat(),
                'timepoints': [],
                'config': {
                    'xdim': self.spin_xdim.value(),
                    'ydim': self.spin_ydim.value(),
                    'n_clusters': self.result.get('n_clusters', 10),
                }
            }
            
            # Ajouter le timepoint actuel
            current_timepoint = {
                'id': f"T{len(patient_data['timepoints'])}",
                'date': datetime.now().strftime("%Y-%m-%d %H:%M"),
                'n_clusters': self.result.get('n_clusters', 0),
                'total_cells': self.result.get('cell_data').shape[0] if self.result.get('cell_data') is not None else 0,
                'files_analyzed': self.pathological_files + self.healthy_files,
                'metacluster_distribution': {}
            }
            
            # Calculer la distribution des metaclusters
            if self.result.get('cell_data') is not None:
                cell_data = self.result['cell_data']
                if 'metaclustering' in cell_data.obs:
                    mc = cell_data.obs['metaclustering'].values
                    total = len(mc)
                    for i in range(self.result.get('n_clusters', 10)):
                        pct = (mc == i).sum() / total * 100
                        current_timepoint['metacluster_distribution'][str(i)] = round(pct, 2)
            
            patient_data['timepoints'].append(current_timepoint)
            
            # Charger les donnees existantes si le patient existe deja
            data_folder = self._get_patient_data_folder()
            patient_file = data_folder / f"{patient_id}.json"
            
            if patient_file.exists():
                try:
                    with open(patient_file, 'r', encoding='utf-8') as f:
                        existing_data = json.load(f)
                    # Fusionner les timepoints
                    patient_data['timepoints'] = existing_data.get('timepoints', []) + patient_data['timepoints']
                    patient_data['created_at'] = existing_data.get('created_at', patient_data['created_at'])
                except Exception:
                    pass
            
            # Sauvegarder
            with open(patient_file, 'w', encoding='utf-8') as f:
                json.dump(patient_data, f, indent=2, ensure_ascii=False)
            
            self._log(f"Patient {patient_id} sauvegarde: {patient_file}")
            self._update_patient_timeline(patient_data)
            
            QMessageBox.information(self, "Succes", f"Patient {patient_id} sauvegarde avec succes!\n{len(patient_data['timepoints'])} timepoint(s)")
            
        except Exception as e:
            QMessageBox.critical(self, "Erreur", f"Erreur lors de la sauvegarde:\n{str(e)}")
            self._log(f"Erreur sauvegarde patient: {str(e)}")
    
    def _load_patient_data(self):
        """Charge les donnees d'un patient existant."""
        data_folder = self._get_patient_data_folder()
        
        # Lister les patients disponibles
        patient_files = list(data_folder.glob("*.json"))
        
        if not patient_files:
            QMessageBox.information(self, "Info", "Aucun patient enregistre.\nLes patients seront stockes dans:\n" + str(data_folder))
            return
        
        # Permettre de selectionner un fichier
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Charger un Patient", str(data_folder),
            "Fichiers Patient (*.json)"
        )
        
        if not file_path:
            return
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                patient_data = json.load(f)
            
            patient_id = patient_data.get('patient_id', Path(file_path).stem)
            
            if self.edit_patient_id is not None:
                self.edit_patient_id.setText(patient_id)
            
            self._update_patient_timeline(patient_data)
            self._plot_patient_evolution(patient_data)
            
            n_timepoints = len(patient_data.get('timepoints', []))
            self._log(f"Patient {patient_id} charge: {n_timepoints} timepoint(s)")
            
            QMessageBox.information(self, "Patient Charge", f"Patient: {patient_id}\nTimepoints: {n_timepoints}")
            
            # Activer le bouton ajouter timepoint
            if self.btn_add_timepoint is not None:
                self.btn_add_timepoint.setEnabled(True)
            
        except Exception as e:
            QMessageBox.critical(self, "Erreur", f"Erreur lors du chargement:\n{str(e)}")
            self._log(f"Erreur chargement patient: {str(e)}")
    
    def _add_patient_timepoint(self):
        """Ajoute le resultat actuel comme nouveau timepoint pour le patient charge."""
        if not self.result:
            QMessageBox.warning(self, "Attention", "Lancez d'abord une analyse FlowSOM.")
            return
            
        if self.edit_patient_id is None:
            return
            
        patient_id = self.edit_patient_id.text().strip()
        if not patient_id:
            QMessageBox.warning(self, "Attention", "Aucun patient charge. Chargez d'abord un patient.")
            return
        
        # Sauvegarder le nouveau timepoint
        self._save_patient_data()
    
    def _update_patient_timeline(self, patient_data: Dict[str, Any]):
        """Met a jour la timeline avec les donnees patient."""
        if self.patient_timeline is None:
            return
        
        timepoints = patient_data.get('timepoints', [])
        self.patient_timeline.setRowCount(len(timepoints))
        
        for i, tp in enumerate(timepoints):
            self.patient_timeline.setItem(i, 0, QTableWidgetItem(tp.get('id', f'T{i}')))
            self.patient_timeline.setItem(i, 1, QTableWidgetItem(tp.get('date', 'N/A')))
            
            files = tp.get('files_analyzed', [])
            file_str = Path(files[0]).name if files else 'N/A'
            self.patient_timeline.setItem(i, 2, QTableWidgetItem(file_str))
            
            self.patient_timeline.setItem(i, 3, QTableWidgetItem(str(tp.get('n_clusters', 'N/A'))))
            self.patient_timeline.setItem(i, 4, QTableWidgetItem(tp.get('notes', '')))
    
    def _plot_patient_evolution(self, patient_data: Dict[str, Any]):
        """Trace l'evolution temporelle du patient."""
        if self.patient_canvas is None:
            return
        
        timepoints = patient_data.get('timepoints', [])
        if len(timepoints) < 1:
            return
        
        self.patient_canvas.clear_figure()
        ax = self.patient_canvas.fig.add_subplot(111)
        
        # Preparer les donnees pour le plot
        n_clusters = max([len(tp.get('metacluster_distribution', {})) for tp in timepoints])
        
        if n_clusters == 0:
            ax.text(0.5, 0.5, "Pas de donnees de distribution disponibles",
                   transform=ax.transAxes, ha='center', va='center', color='#a6adc8')
            ax.set_facecolor('#1e1e2e')
            self.patient_canvas.fig.patch.set_facecolor('#1e1e2e')
            self.patient_canvas.draw()
            return
        
        # Couleurs pour chaque cluster
        colors = plt.cm.tab20(np.linspace(0, 1, n_clusters))
        
        x = range(len(timepoints))
        labels = [tp.get('id', f'T{i}') for i, tp in enumerate(timepoints)]
        
        # Stacked bar chart
        bottom = np.zeros(len(timepoints))
        
        for cluster_idx in range(n_clusters):
            values = []
            for tp in timepoints:
                dist = tp.get('metacluster_distribution', {})
                values.append(dist.get(str(cluster_idx), 0))
            
            ax.bar(x, values, bottom=bottom, label=f'MC{cluster_idx}', 
                  color=colors[cluster_idx], edgecolor='white', linewidth=0.3)
            bottom += np.array(values)
        
        ax.set_xticks(x)
        ax.set_xticklabels(labels, color='#cdd6f4')
        ax.set_ylabel('Distribution (%)', color='#cdd6f4')
        ax.set_xlabel('Timepoint', color='#cdd6f4')
        ax.set_title(f"Evolution Patient: {patient_data.get('patient_id', 'N/A')}", 
                    color='#cdd6f4', fontsize=12, fontweight='bold')
        ax.tick_params(colors='#cdd6f4')
        ax.set_facecolor('#1e1e2e')
        ax.spines['bottom'].set_color('#45475a')
        ax.spines['left'].set_color('#45475a')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Legende
        if n_clusters <= 10:
            ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=8,
                     facecolor='#313244', edgecolor='#45475a')
        
        self.patient_canvas.fig.patch.set_facecolor('#1e1e2e')
        self.patient_canvas.fig.tight_layout()
        self.patient_canvas.draw()
        
    def _run_analysis(self):
        if not FLOWSOM_AVAILABLE:
            QMessageBox.critical(self, "Erreur", "FlowSOM n'est pas installé.\npip install flowsom")
            return
            
        if not self.healthy_files and not self.pathological_files:
            QMessageBox.warning(self, "Attention", "Chargez au moins un dossier FCS.")
            return
            
        params = {
            'xdim': self.spin_xdim.value(),
            'ydim': self.spin_ydim.value(),
            'n_clusters': self.spin_n_clusters.value(),
            'seed': self.spin_seed.value(),
            'auto_cluster': self.chk_auto_cluster.isChecked(),
            'exclude_scatter': self.chk_exclude_scatter.isChecked(),
            'pregate': self.chk_pregate.isChecked(),
            'transform': self.combo_transform.currentText(),
            'cofactor': self.spin_cofactor.value()
        }
        
        self._log("=" * 50)
        self._log("Demarrage de l'analyse FlowSOM")
        self._log(f"   Transformation: {params['transform']}")
        self._log(f"   Pre-gating: {'Oui' if params['pregate'] else 'Non'}")
        self._log(f"   Grille: {params['xdim']}x{params['ydim']}")
        self._log(f"   Metaclusters: {params['n_clusters']}")
        
        self.btn_run.setEnabled(False)
        self.btn_cancel.setEnabled(True)
        self.btn_healthy.setEnabled(False)
        self.btn_pathological.setEnabled(False)
        
        self.worker = FlowSOMWorker(self.healthy_files, self.pathological_files, params)
        self.worker.progress.connect(self._on_progress)
        self.worker.status.connect(self._on_status)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.start()
        
    def _cancel_analysis(self):
        if self.worker and self.worker.isRunning():
            self.worker.cancel()
            self.worker.terminate()
            self.worker.wait()
            self._log("Analyse annulee")
            self._reset_ui_state()
            
    def _on_progress(self, value):
        self.progress_bar.setValue(value)
        
    def _on_status(self, message):
        self.lbl_status.setText(message)
        self.statusBar().showMessage(message)
        self._log(message)
        
    def _on_finished(self, result):
        self.result = result
        self._log("Analyse terminee!")
        
        # Mettre a jour l'interface immediatement
        self.lbl_status.setText("Preparation des resultats...")
        QApplication.processEvents()
        
        self.combo_marker.clear()
        self.combo_marker.setEnabled(True)
        
        var_names = result.get('var_names', [])
        cols_to_use = result.get('cols_to_use', [])
        used_markers = [var_names[i] for i in cols_to_use if i < len(var_names)]
        self.combo_marker.addItems(used_markers)
        
        self.btn_export_fcs.setEnabled(True)
        self.btn_export_csv.setEnabled(True)
        self.btn_export_fig.setEnabled(True)
        self.btn_export_pdf.setEnabled(True)
        self.btn_save_patient.setEnabled(True)
        
        # Remplir la liste des clusters (signaux bloques, pas de selection auto)
        self._populate_cluster_list()
        
        # Activer le bouton ajouter timepoint
        if self.btn_add_timepoint is not None:
            self.btn_add_timepoint.setEnabled(True)
        
        # Reset UI state d'abord pour reactiver les boutons
        self._reset_ui_state()
        
        # Afficher un message de succes immediatement
        self.lbl_status.setText("Analyse terminee! Generation des statistiques...")
        self.lbl_status.setStyleSheet("color: #a6e3a1; padding: 8px; font-size: 10pt; font-weight: 600;")
        QApplication.processEvents()
        
        # Mettre a jour les statistiques (leger)
        self._update_statistics()
        QApplication.processEvents()
        
        # Mettre a jour la visualisation (peut etre lent selon le type)
        self.lbl_status.setText("Generation de la visualisation...")
        QApplication.processEvents()
        self._update_visualization()
        
        self.lbl_status.setText("Analyse terminee avec succes!")
        
    def _on_error(self, error_message):
        self._log(f"ERREUR: {error_message}")
        QMessageBox.critical(self, "Erreur", error_message)
        self._reset_ui_state()
        self.lbl_status.setText("Erreur lors de l'analyse")
        self.lbl_status.setStyleSheet("color: #f38ba8; padding: 8px; font-size: 10pt; font-weight: 600;")
        
    def _reset_ui_state(self):
        self.btn_run.setEnabled(True)
        self.btn_cancel.setEnabled(False)
        self.btn_healthy.setEnabled(True)
        self.btn_pathological.setEnabled(True)
        self.progress_bar.setValue(0)
        
    def _update_statistics(self):
        if not self.result:
            return
            
        cell_data = self.result['cell_data']
        n_clusters = self.result['n_clusters']
        
        metaclusters = cell_data.obs['metaclustering']
        conditions = cell_data.obs['condition']
        
        self.table_stats.setRowCount(n_clusters)
        
        for i in range(n_clusters):
            mask_cluster = metaclusters == i
            
            mask_healthy = (conditions == 'Sain') & mask_cluster
            total_healthy = (conditions == 'Sain').sum()
            pct_h = (mask_healthy.sum() / total_healthy * 100) if total_healthy > 0 else 0
            
            mask_patho = (conditions == 'Pathologique') & mask_cluster
            total_patho = (conditions == 'Pathologique').sum()
            pct_p = (mask_patho.sum() / total_patho * 100) if total_patho > 0 else 0
            
            diff = pct_p - pct_h
            
            self.table_stats.setItem(i, 0, QTableWidgetItem(f"Cluster {i}"))
            self.table_stats.setItem(i, 1, QTableWidgetItem(f"{pct_h:.2f}%"))
            self.table_stats.setItem(i, 2, QTableWidgetItem(f"{pct_p:.2f}%"))
            
            diff_item = QTableWidgetItem(f"{diff:+.2f}%")
            if diff > 5:
                diff_item.setForeground(QColor("#f38ba8"))
            elif diff < -5:
                diff_item.setForeground(QColor("#a6e3a1"))
            self.table_stats.setItem(i, 3, diff_item)
            
    def _update_visualization(self):
        if not self.result:
            return
            
        viz_type = self.combo_viz_type.currentText()
        self.combo_marker.setEnabled("Marker" in viz_type)
        self.canvas.clear_figure()
        
        try:
            if "Star Chart" in viz_type:
                self._plot_star_chart()
            elif "Arbre MST" in viz_type:
                self._plot_mst_tree_matplotlib()
            elif "Grille SOM" in viz_type:
                self._plot_som_grid_matplotlib()
            elif "Grid View" in viz_type:
                self._plot_grid_view()
            elif "Cluster Numbers" in viz_type:
                self._plot_cluster_numbers()
            elif "Heatmap" in viz_type:
                self._plot_heatmap()
            elif "Distribution" in viz_type:
                self._plot_condition_distribution()
            elif "Marker" in viz_type:
                self._plot_marker_expression()
            elif "UMAP" in viz_type:
                self._plot_umap()
            elif "t-SNE" in viz_type:
                self._plot_tsne()
        except Exception as e:
            self._log(f"Erreur visualisation: {str(e)}")
            ax = self.canvas.fig.add_subplot(111)
            ax.text(0.5, 0.5, f"Erreur: {str(e)}", transform=ax.transAxes, 
                   ha='center', va='center', color='#f38ba8', fontsize=12)
            ax.set_facecolor('#1e1e2e')
            
        self.canvas.draw()
        
    def _plot_star_chart(self):
        """Affiche le Star Chart avec l'API FlowSOM."""
        fsom = self.result['fsom']
        
        try:
            fig = fs.pl.plot_stars(
                fsom,
                background_values=fsom.get_cluster_data().obs.metaclustering,
                view="MST"
            )
            self._embed_flowsom_figure(fig, "FlowSOM - Star Chart (MST)")
        except Exception as e:
            self._fallback_visualization(f"Star Chart: {str(e)}")
            
    def _plot_grid_view(self):
        """Affiche la Grid View avec l'API FlowSOM."""
        fsom = self.result['fsom']
        
        try:
            fig = fs.pl.plot_stars(
                fsom,
                background_values=fsom.get_cluster_data().obs.metaclustering,
                view="grid",
                equal_node_size=True,
                equal_background_size=True
            )
            self._embed_flowsom_figure(fig, "FlowSOM - Grid View")
        except Exception as e:
            self._fallback_visualization(f"Grid View: {str(e)}")
    
    def _plot_mst_tree_matplotlib(self):
        """
        Affiche l'arbre MST (Minimum Spanning Tree) en matplotlib.
        Utilise les coordonnees layout (xNodes, yNodes) du FlowSOM.
        """
        fsom = self.result['fsom']
        cluster_data = fsom.get_cluster_data()
        cell_data = fsom.get_cell_data()
        
        try:
            # Recuperer les coordonnees du layout MST
            layout = cluster_data.obsm.get('layout', None)
            if layout is None:
                self._fallback_visualization("Coordonnees MST non disponibles")
                return
            
            # Recuperer le clustering et metaclustering
            clustering = cell_data.obs['clustering'].values
            metaclustering = cluster_data.obs['metaclustering'].values
            
            # Calculer la taille de chaque node (nombre de cellules)
            n_nodes = len(cluster_data)
            node_sizes = np.zeros(n_nodes)
            for i in range(n_nodes):
                node_sizes[i] = np.sum(clustering == i)
            
            # Normaliser les tailles pour l'affichage
            max_size = node_sizes.max() if node_sizes.max() > 0 else 1
            sizes = 50 + (node_sizes / max_size) * 500  # Entre 50 et 550
            
            # Palette de couleurs pour les metaclusters
            n_meta = len(np.unique(metaclustering))
            cmap = plt.cm.tab20 if n_meta <= 20 else plt.cm.turbo
            colors = [cmap(int(m) / max(n_meta - 1, 1)) for m in metaclustering]
            
            ax = self.canvas.fig.add_subplot(111)
            
            # Dessiner les aretes du MST si disponibles
            try:
                # Essayer de recuperer les aretes du MST depuis l'objet FlowSOM
                if hasattr(fsom, 'mudata') and 'cluster_data' in fsom.mudata.mod:
                    cd = fsom.mudata.mod['cluster_data']
                    if 'mst' in cd.uns:
                        mst = cd.uns['mst']
                        # Dessiner les aretes
                        for edge in mst:
                            i, j = int(edge[0]), int(edge[1])
                            ax.plot([layout[i, 0], layout[j, 0]], 
                                   [layout[i, 1], layout[j, 1]],
                                   color='#6c7086', linewidth=1, alpha=0.6, zorder=1)
            except:
                # Si pas d'aretes MST, connecter les nodes proches
                from scipy.spatial.distance import pdist, squareform
                distances = squareform(pdist(layout))
                threshold = np.percentile(distances[distances > 0], 15)
                for i in range(n_nodes):
                    for j in range(i + 1, n_nodes):
                        if distances[i, j] < threshold:
                            ax.plot([layout[i, 0], layout[j, 0]], 
                                   [layout[i, 1], layout[j, 1]],
                                   color='#6c7086', linewidth=0.8, alpha=0.4, zorder=1)
            
            # Scatter plot des nodes
            scatter = ax.scatter(layout[:, 0], layout[:, 1], 
                               s=sizes, c=colors, edgecolors='white', 
                               linewidths=1.5, alpha=0.9, zorder=2)
            
            # Annoter avec les numeros de metaclusters
            for i in range(n_nodes):
                if node_sizes[i] > 0:  # N'afficher que les nodes non vides
                    ax.annotate(str(int(metaclustering[i])), 
                               (layout[i, 0], layout[i, 1]),
                               ha='center', va='center', fontsize=7,
                               color='#11111b', fontweight='bold', zorder=3)
            
            # Style
            ax.set_xlabel('xNodes', color='#cdd6f4', fontsize=12, fontweight='bold')
            ax.set_ylabel('yNodes', color='#cdd6f4', fontsize=12, fontweight='bold')
            ax.set_title(f'Arbre MST - {n_nodes} nodes, {n_meta} metaclusters', 
                        fontsize=14, color='#cdd6f4', pad=15, fontweight='bold')
            ax.tick_params(colors='#cdd6f4', labelsize=10)
            ax.set_facecolor('#1e1e2e')
            ax.spines['bottom'].set_color('#45475a')
            ax.spines['left'].set_color('#45475a')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(True, alpha=0.15, color='#6c7086', linestyle='--')
            
            # Legende des metaclusters
            from matplotlib.patches import Patch
            if n_meta <= 15:
                legend_handles = [Patch(facecolor=cmap(i / max(n_meta - 1, 1)), 
                                       edgecolor='white', label=f'MC {i}') 
                                 for i in range(n_meta)]
                ax.legend(handles=legend_handles, loc='upper right', 
                         fontsize=7, facecolor='#313244', labelcolor='#cdd6f4',
                         edgecolor='#45475a', framealpha=0.9, ncol=2)
            
            self.canvas.fig.patch.set_facecolor('#1e1e2e')
            self.canvas.fig.tight_layout(pad=1.5)
            
        except Exception as e:
            self._fallback_visualization(f"Arbre MST: {str(e)}")
    
    def _plot_som_grid_matplotlib(self):
        """
        Affiche la grille SOM en matplotlib.
        Utilise les coordonnees grid (xGrid, yGrid) du FlowSOM.
        """
        fsom = self.result['fsom']
        cluster_data = fsom.get_cluster_data()
        cell_data = fsom.get_cell_data()
        
        try:
            # Recuperer les coordonnees de la grille
            grid = cluster_data.obsm.get('grid', None)
            if grid is None:
                self._fallback_visualization("Coordonnees grille non disponibles")
                return
            
            # Recuperer le clustering et metaclustering
            clustering = cell_data.obs['clustering'].values
            metaclustering = cluster_data.obs['metaclustering'].values
            
            # Calculer la taille de chaque node
            n_nodes = len(cluster_data)
            node_sizes = np.zeros(n_nodes)
            for i in range(n_nodes):
                node_sizes[i] = np.sum(clustering == i)
            
            # Palette de couleurs pour les metaclusters
            n_meta = len(np.unique(metaclustering))
            cmap = plt.cm.tab20 if n_meta <= 20 else plt.cm.turbo
            colors = [cmap(int(m) / max(n_meta - 1, 1)) for m in metaclustering]
            
            # Normaliser les tailles - plus grand pour la grille
            max_size = node_sizes.max() if node_sizes.max() > 0 else 1
            sizes = 100 + (node_sizes / max_size) * 800  # Plus grands pour la grille
            
            ax = self.canvas.fig.add_subplot(111)
            
            # Scatter plot avec taille proportionnelle
            scatter = ax.scatter(grid[:, 0], grid[:, 1], 
                               s=sizes, c=colors, edgecolors='white', 
                               linewidths=2, alpha=0.85, marker='o')
            
            # Ajouter les numeros de cluster et metacluster
            for i in range(n_nodes):
                if node_sizes[i] > 0:
                    # Afficher le metacluster en gras
                    ax.annotate(str(int(metaclustering[i])), 
                               (grid[i, 0], grid[i, 1]),
                               ha='center', va='center', fontsize=9,
                               color='#11111b', fontweight='bold')
            
            # Determiner les dimensions de la grille
            xdim = int(grid[:, 0].max()) + 1
            ydim = int(grid[:, 1].max()) + 1
            
            # Style
            ax.set_xlabel('xGrid', color='#cdd6f4', fontsize=12, fontweight='bold')
            ax.set_ylabel('yGrid', color='#cdd6f4', fontsize=12, fontweight='bold')
            ax.set_title(f'Grille SOM {xdim}x{ydim} - {n_nodes} nodes, {n_meta} metaclusters', 
                        fontsize=14, color='#cdd6f4', pad=15, fontweight='bold')
            ax.tick_params(colors='#cdd6f4', labelsize=10)
            ax.set_facecolor('#1e1e2e')
            ax.spines['bottom'].set_color('#45475a')
            ax.spines['left'].set_color('#45475a')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            # Grille de fond
            ax.set_xticks(np.arange(xdim))
            ax.set_yticks(np.arange(ydim))
            ax.grid(True, alpha=0.2, color='#6c7086', linestyle='-')
            
            # Ajuster les limites pour centrer
            ax.set_xlim(-0.5, xdim - 0.5)
            ax.set_ylim(-0.5, ydim - 0.5)
            ax.set_aspect('equal')
            
            # Legende des metaclusters
            from matplotlib.patches import Patch
            if n_meta <= 15:
                legend_handles = [Patch(facecolor=cmap(i / max(n_meta - 1, 1)), 
                                       edgecolor='white', label=f'MC {i}') 
                                 for i in range(n_meta)]
                ax.legend(handles=legend_handles, loc='upper right', 
                         fontsize=7, facecolor='#313244', labelcolor='#cdd6f4',
                         edgecolor='#45475a', framealpha=0.9, ncol=2)
            
            # Ajouter une barre de couleur pour les tailles
            # (optionnel - commenté pour garder simple)
            # from mpl_toolkits.axes_grid1 import make_axes_locatable
            
            self.canvas.fig.patch.set_facecolor('#1e1e2e')
            self.canvas.fig.tight_layout(pad=1.5)
            
        except Exception as e:
            self._fallback_visualization(f"Grille SOM: {str(e)}")
            
    def _plot_cluster_numbers(self):
        """Affiche les numéros de clusters."""
        fsom = self.result['fsom']
        
        try:
            fig = fs.pl.plot_numbers(fsom, level="metaclusters", text_size=8)
            self._embed_flowsom_figure(fig, "FlowSOM - Numéros de Clusters")
        except Exception as e:
            self._fallback_visualization(f"Cluster Numbers: {str(e)}")
            
    def _plot_marker_expression(self):
        """Affiche l'expression d'un marqueur."""
        marker = self.combo_marker.currentText()
        if not marker:
            return
            
        fsom = self.result['fsom']
        
        try:
            fig = fs.pl.plot_marker(fsom, marker=np.array([marker]))
            self._embed_flowsom_figure(fig, f"Expression: {marker}")
        except Exception as e:
            self._fallback_visualization(f"Marker Expression: {str(e)}")
            
    def _embed_flowsom_figure(self, fig, title):
        """Intègre une figure FlowSOM dans le canvas."""
        if fig is None:
            return
            
        try:
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                tmp_path = tmp.name
            
            fig.savefig(tmp_path, dpi=150, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            plt.close(fig)
            
            img = plt.imread(tmp_path)
            ax = self.canvas.fig.add_subplot(111)
            ax.imshow(img)
            ax.axis('off')
            ax.set_title(title, fontsize=14, color='#cdd6f4', pad=15, fontweight='bold')
            self.canvas.fig.patch.set_facecolor('#1e1e2e')
            
            try:
                os.unlink(tmp_path)
            except:
                pass
        except Exception as e:
            self._fallback_visualization(str(e))
            
    def _fallback_visualization(self, error_msg):
        """Affiche un message d'erreur."""
        ax = self.canvas.fig.add_subplot(111)
        ax.text(0.5, 0.5, f"{error_msg}\n\nEssayez une autre visualisation.", 
               transform=ax.transAxes, ha='center', va='center', 
               color='#fab387', fontsize=11, wrap=True)
        ax.set_facecolor('#1e1e2e')
        self.canvas.fig.patch.set_facecolor('#1e1e2e')
            
    def _plot_heatmap(self):
        """Affiche la heatmap des expressions."""
        fsom = self.result['fsom']
        cell_data = fsom.get_cell_data()
        
        n_clusters = self.result['n_clusters']
        cols_to_use = self.result['cols_to_use']
        var_names = self.result['var_names']
        used_markers = [var_names[i] for i in cols_to_use if i < len(var_names)]
        
        X = cell_data.X
        if hasattr(X, 'toarray'):
            X = X.toarray()
        
        metaclustering = cell_data.obs['metaclustering'].values
        
        mfi_matrix = np.zeros((n_clusters, len(cols_to_use)))
        for i in range(n_clusters):
            mask = metaclustering == i
            if mask.sum() > 0:
                mfi_matrix[i, :] = np.nanmean(X[mask][:, cols_to_use], axis=0)
        
        mfi_normalized = (mfi_matrix - np.nanmean(mfi_matrix, axis=0)) / (np.nanstd(mfi_matrix, axis=0) + 1e-10)
        
        ax = self.canvas.fig.add_subplot(111)
        im = ax.imshow(mfi_normalized.T, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
        
        ax.set_yticks(range(len(used_markers)))
        ax.set_yticklabels(used_markers, fontsize=9, color='#cdd6f4')
        ax.set_xticks(range(n_clusters))
        ax.set_xticklabels([f"MC{i}" for i in range(n_clusters)], fontsize=9, color='#cdd6f4')
        
        ax.set_title("Heatmap - Expression par Métacluster (Z-score)", 
                    fontsize=14, color='#cdd6f4', pad=15, fontweight='bold')
        ax.set_xlabel("Métacluster", color='#cdd6f4', fontsize=11)
        ax.set_ylabel("Marqueur", color='#cdd6f4', fontsize=11)
        ax.tick_params(colors='#cdd6f4')
        ax.set_facecolor('#1e1e2e')
        
        cbar = self.canvas.fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.ax.yaxis.set_tick_params(color='#cdd6f4')
        cbar.outline.set_edgecolor('#45475a')
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='#cdd6f4')
        
    def _plot_condition_distribution(self):
        """Affiche la distribution par condition."""
        cell_data = self.result['cell_data']
        n_clusters = self.result['n_clusters']
        
        metaclustering = cell_data.obs['metaclustering'].values
        conditions = cell_data.obs['condition'].values
        
        healthy_pcts = []
        patho_pcts = []
        
        for i in range(n_clusters):
            mask_cluster = metaclustering == i
            
            mask_healthy = (conditions == 'Sain') & mask_cluster
            total_healthy = (conditions == 'Sain').sum()
            healthy_pcts.append((mask_healthy.sum() / total_healthy * 100) if total_healthy > 0 else 0)
            
            mask_patho = (conditions == 'Pathologique') & mask_cluster
            total_patho = (conditions == 'Pathologique').sum()
            patho_pcts.append((mask_patho.sum() / total_patho * 100) if total_patho > 0 else 0)
        
        ax = self.canvas.fig.add_subplot(111)
        x = np.arange(n_clusters)
        width = 0.35
        
        bars1 = ax.bar(x - width/2, healthy_pcts, width, label='NBM', color='#a6e3a1', edgecolor='white', linewidth=0.5)
        bars2 = ax.bar(x + width/2, patho_pcts, width, label='Patient', color='#f38ba8', edgecolor='white', linewidth=0.5)
        
        ax.set_xlabel('Métacluster', color='#cdd6f4', fontsize=11)
        ax.set_ylabel('Pourcentage (%)', color='#cdd6f4', fontsize=11)
        ax.set_title('Distribution des Métaclusters', fontsize=14, color='#cdd6f4', pad=15, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels([f'MC{i}' for i in range(n_clusters)], color='#cdd6f4', fontsize=9)
        ax.tick_params(colors='#cdd6f4')
        ax.legend(facecolor='#313244', labelcolor='#cdd6f4', edgecolor='#45475a')
        ax.set_facecolor('#1e1e2e')
        ax.spines['bottom'].set_color('#45475a')
        ax.spines['left'].set_color('#45475a')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.2, color='#45475a')
        
    def _export_fcs(self):
        if not self.result:
            return
            
        filepath, _ = QFileDialog.getSaveFileName(self, "Exporter FCS", "flowsom_results.fcs", "FCS (*.fcs)")
        
        if filepath:
            try:
                cell_data = self.result['cell_data']
                X = cell_data.X
                if hasattr(X, 'toarray'):
                    X = X.toarray()
                
                self._log(f"Export FCS: {X.shape[0]} cellules, {X.shape[1]} parametres originaux")
                
                # Recuperer les clusters et metaclusters
                clustering = cell_data.obs['clustering'].values
                metaclustering = cell_data.obs['metaclustering'].values
                
                self._log(f"   Clusters uniques: {len(np.unique(clustering))}")
                self._log(f"   Metaclusters uniques: {len(np.unique(metaclustering))}")
                
                # Convertir en float32 et reshape
                clustering_col = clustering.astype(np.float32).reshape(-1, 1)
                metaclustering_col = metaclustering.astype(np.float32).reshape(-1, 1)
                
                # Preparer les donnees de base
                X_float = X.astype(np.float32)
                export_data = np.hstack([X_float, clustering_col, metaclustering_col])
                channel_names = list(cell_data.var_names) + ['FlowSOM_cluster', 'FlowSOM_metacluster']
                
                self._log(f"   Donnees preparees: {export_data.shape}")
                
                # Ajouter condition et file_origin comme colonnes numeriques
                # Encoder les conditions: Sain=0, Pathologique=1
                conditions = cell_data.obs['condition'].values
                condition_encoded = np.array([0.0 if c == 'Sain' else 1.0 for c in conditions], dtype=np.float32).reshape(-1, 1)
                export_data = np.hstack([export_data, condition_encoded])
                channel_names.append('Condition')
                
                # =====================================================================
                # AJOUTER LES COORDONNEES SOM (comme dans Kaluza)
                # =====================================================================
                fsom = self.result.get('fsom')
                if fsom is not None:
                    try:
                        cluster_data = fsom.get_cluster_data()
                        n_cells = X.shape[0]
                        
                        # Recuperer les coordonnees de grille (xGrid, yGrid) pour chaque node
                        grid_coords = cluster_data.obsm.get('grid', None)  # shape (n_nodes, 2)
                        # Recuperer les coordonnees du layout MST (xNodes, yNodes)
                        layout_coords = cluster_data.obsm.get('layout', None)  # shape (n_nodes, 2)
                        
                        # Calculer la taille de chaque node (nombre de cellules)
                        node_sizes = np.zeros(len(cluster_data), dtype=np.float32)
                        for i in range(len(cluster_data)):
                            node_sizes[i] = np.sum(clustering == i)
                        
                        # Mapper les coordonnees sur chaque cellule via son cluster assignment
                        if grid_coords is not None:
                            xGrid = np.array([grid_coords[int(c), 0] for c in clustering], dtype=np.float32).reshape(-1, 1)
                            yGrid = np.array([grid_coords[int(c), 1] for c in clustering], dtype=np.float32).reshape(-1, 1)
                            export_data = np.hstack([export_data, xGrid, yGrid])
                            channel_names.extend(['xGrid', 'yGrid'])
                            self._log(f"   Coordonnees grille ajoutees (xGrid, yGrid)")
                        
                        if layout_coords is not None:
                            xNodes = np.array([layout_coords[int(c), 0] for c in clustering], dtype=np.float32).reshape(-1, 1)
                            yNodes = np.array([layout_coords[int(c), 1] for c in clustering], dtype=np.float32).reshape(-1, 1)
                            export_data = np.hstack([export_data, xNodes, yNodes])
                            channel_names.extend(['xNodes', 'yNodes'])
                            self._log(f"   Coordonnees MST ajoutees (xNodes, yNodes)")
                        
                        # Ajouter la taille du node pour chaque cellule
                        size_col = np.array([node_sizes[int(c)] for c in clustering], dtype=np.float32).reshape(-1, 1)
                        export_data = np.hstack([export_data, size_col])
                        channel_names.append('size')
                        self._log(f"   Taille des nodes ajoutee (size)")
                        
                    except Exception as e:
                        self._log(f"   ATTENTION: Impossible d'ajouter les coordonnees SOM: {e}")
                
                # =====================================================================
                # AJOUTER UMAP si disponible
                # =====================================================================
                # Utilise -999 comme valeur manquante (facilement filtrable)
                MISSING_VALUE = -999.0
                umap_cells = 0
                if 'X_umap' in cell_data.obsm:
                    umap_coords = cell_data.obsm['X_umap'].astype(np.float32).copy()
                    umap_cells = np.sum(~np.isnan(umap_coords[:, 0]))
                    # Remplacer NaN par valeur manquante avant export
                    umap_coords = np.nan_to_num(umap_coords, nan=MISSING_VALUE)
                    export_data = np.hstack([export_data, umap_coords])
                    channel_names.extend(['UMAP1', 'UMAP2'])
                    self._log(f"   UMAP coordonnees ajoutees ({umap_cells:,} / {X.shape[0]:,} cellules valides)")
                
                # Ajouter t-SNE si disponible
                tsne_cells = 0
                if 'X_tsne' in cell_data.obsm:
                    tsne_coords = cell_data.obsm['X_tsne'].astype(np.float32).copy()
                    tsne_cells = np.sum(~np.isnan(tsne_coords[:, 0]))
                    # Remplacer NaN par valeur manquante avant export
                    tsne_coords = np.nan_to_num(tsne_coords, nan=MISSING_VALUE)
                    export_data = np.hstack([export_data, tsne_coords])
                    channel_names.extend(['tSNE1', 'tSNE2'])
                    self._log(f"   t-SNE coordonnees ajoutees ({tsne_cells:,} / {X.shape[0]:,} cellules valides)")
                
                self._log(f"   Export final: {export_data.shape[0]} events, {export_data.shape[1]} canaux")
                self._log(f"   Canaux: {channel_names}")
                
                # Note: Les coordonnees UMAP/t-SNE manquantes sont deja remplacees par -999
                # Verifier les autres NaN ou Inf
                nan_count = np.sum(np.isnan(export_data))
                if nan_count > 0:
                    self._log(f"   ATTENTION: {nan_count} NaN restants, remplacement par 0")
                    export_data = np.nan_to_num(export_data, nan=0.0)
                inf_count = np.sum(np.isinf(export_data))
                if inf_count > 0:
                    self._log(f"   ATTENTION: {inf_count} Inf detectes, remplacement par max")
                    export_data = np.nan_to_num(export_data, posinf=1e6, neginf=-1e6)
                
                # Utiliser la methode d'ecriture FCS
                self._write_fcs_with_flowio(filepath, channel_names, export_data)
                
                self._log(f"FCS exporte avec succes: {filepath}")
                
                # Construire le message de succes detaille
                msg_details = [
                    f"Fichier exporte:\n{filepath}\n",
                    f"Cellules: {export_data.shape[0]:,}",
                    f"Canaux: {len(channel_names)}",
                    "",
                    "Inclus:",
                    "  - Tous les marqueurs originaux",
                    "  - FlowSOM_cluster (nodes SOM)",
                    "  - FlowSOM_metacluster (metaclusters)",
                    "  - Condition (0=NBM, 1=Patient)",
                    "  - xGrid, yGrid (coordonnees grille SOM)",
                    "  - xNodes, yNodes (coordonnees MST)",
                    "  - size (taille du node)"
                ]
                if 'X_umap' in cell_data.obsm:
                    msg_details.append(f"  - UMAP1, UMAP2 ({umap_cells:,} valides)")
                if 'X_tsne' in cell_data.obsm:
                    msg_details.append(f"  - tSNE1, tSNE2 ({tsne_cells:,} valides)")
                
                QMessageBox.information(self, "Export FCS Reussi", "\n".join(msg_details))
            except Exception as e:
                import traceback
                self._log(f"Erreur export: {str(e)}")
                self._log(traceback.format_exc())
                QMessageBox.critical(self, "Erreur", str(e))
                
    def _export_csv(self):
        if not self.result:
            return
            
        filepath, _ = QFileDialog.getSaveFileName(self, "Exporter CSV", "flowsom_results.csv", "CSV (*.csv)")
        
        if filepath:
            try:
                cell_data = self.result['cell_data']
                X = cell_data.X
                if hasattr(X, 'toarray'):
                    X = X.toarray()
                
                self._log(f"Export CSV: {X.shape[0]} cellules")
                    
                df = pd.DataFrame(X, columns=list(cell_data.var_names))
                df['condition'] = cell_data.obs['condition'].values
                df['file_origin'] = cell_data.obs['file_origin'].values
                df['FlowSOM_cluster'] = cell_data.obs['clustering'].values
                df['FlowSOM_metacluster'] = cell_data.obs['metaclustering'].values
                
                # Ajouter UMAP si disponible
                if 'X_umap' in cell_data.obsm:
                    umap_coords = cell_data.obsm['X_umap']
                    df['UMAP1'] = umap_coords[:, 0]
                    df['UMAP2'] = umap_coords[:, 1]
                    self._log("   UMAP coordonnees ajoutees")
                
                # Ajouter t-SNE si disponible
                if 'X_tsne' in cell_data.obsm:
                    tsne_coords = cell_data.obsm['X_tsne']
                    df['tSNE1'] = tsne_coords[:, 0]
                    df['tSNE2'] = tsne_coords[:, 1]
                    self._log("   t-SNE coordonnees ajoutees")
                
                df.to_csv(filepath, index=False)
                
                self._log(f"CSV exporte: {filepath}")
                self._log(f"   Lignes: {len(df)}, Colonnes: {len(df.columns)}")
                QMessageBox.information(self, "Succes", 
                    f"Fichier exporte:\n{filepath}\n\n"
                    f"Cellules: {len(df):,}\n"
                    f"Colonnes: {len(df.columns)}")
            except Exception as e:
                import traceback
                self._log(f"Erreur export CSV: {str(e)}")
                self._log(traceback.format_exc())
                QMessageBox.critical(self, "Erreur", str(e))
            except Exception as e:
                self._log(f"Erreur export: {str(e)}")
                QMessageBox.critical(self, "Erreur", str(e))
                
    def _export_figure(self):
        if not self.result:
            return
            
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Exporter Figure", "flowsom_figure.png",
            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)"
        )
        
        if filepath:
            try:
                self.canvas.fig.savefig(filepath, dpi=800, bbox_inches='tight', facecolor='#1e1e2e')
                self._log(f"Figure exportee: {filepath}")
                QMessageBox.information(self, "Succes", f"Figure exportee:\n{filepath}")
            except Exception as e:
                self._log(f"Erreur export: {str(e)}")
                QMessageBox.critical(self, "Erreur", str(e))
    
    def _write_fcs_with_flowio(self, filepath: str, channel_names: List[str], data: np.ndarray):
        """
        Ecrit un fichier FCS valide avec toutes les metadonnees requises.
        Utilise le format FCS 3.1 avec les keywords correctement definis.
        """
        n_events, n_channels = data.shape
        self._log(f"Ecriture FCS: {n_events} events, {n_channels} channels")
        
        # S'assurer que les donnees sont en float32 et contiguës en memoire
        data = np.ascontiguousarray(data.astype(np.float32))
        
        try:
            import fcswrite
            # Utiliser fcswrite avec les parametres nommes explicites
            fcswrite.write_fcs(
                filename=filepath, 
                chn_names=channel_names,
                data=data,
                compat_chn_names=True,  # Nettoyer les noms de canaux
                compat_copy=True  # Creer une copie pour eviter les problemes de memoire
            )
            self._log(f"FCS ecrit avec fcswrite: {data.shape[0]} events, {data.shape[1]} channels")
            return
        except ImportError:
            self._log("fcswrite non disponible, utilisation methode manuelle")
        except Exception as e:
            self._log(f"fcswrite echoue: {e}, utilisation methode manuelle")
        
        # Methode manuelle si fcswrite n'est pas disponible
        import struct
        
        n_events, n_channels = data.shape
        data = data.astype(np.float32)
        
        # Construire le TEXT segment avec toutes les metadonnees requises
        # Utiliser un dictionnaire ordonne pour garantir l'ordre
        from collections import OrderedDict
        text_dict = OrderedDict()
        
        # Parametres globaux (obligatoires)
        text_dict['$BEGINANALYSIS'] = '0'
        text_dict['$ENDANALYSIS'] = '0'
        text_dict['$BEGINSTEXT'] = '0'
        text_dict['$ENDSTEXT'] = '0'
        text_dict['$BYTEORD'] = '1,2,3,4'  # Little endian
        text_dict['$DATATYPE'] = 'F'  # Float
        text_dict['$MODE'] = 'L'  # List mode
        text_dict['$NEXTDATA'] = '0'
        text_dict['$PAR'] = str(n_channels)
        text_dict['$TOT'] = str(n_events)
        
        # Ajouter les parametres pour chaque canal
        # Ordre important: N, S, R, B, E pour chaque canal avant de passer au suivant
        for i, name in enumerate(channel_names, 1):
            # Nettoyer le nom du canal
            clean_name = str(name).replace('/', '-').replace('\\', '-').replace('$', '')
            clean_name = clean_name[:64]  # Limiter la longueur
            
            # Calculer le range (max value)
            col_data = data[:, i-1]
            col_data_clean = col_data[~np.isnan(col_data) & ~np.isinf(col_data)]
            if len(col_data_clean) > 0:
                max_val = max(float(np.max(col_data_clean)), 1.0)
                min_val = min(float(np.min(col_data_clean)), 0.0)
            else:
                max_val = 262144.0
                min_val = 0.0
            
            # Range doit etre >= max_val
            range_val = int(np.ceil(max_val)) + 1
            
            # Parametres OBLIGATOIRES dans l'ordre FCS standard
            text_dict[f'$P{i}N'] = clean_name  # Short name (OBLIGATOIRE)
            text_dict[f'$P{i}S'] = clean_name  # Long name (OBLIGATOIRE)
            text_dict[f'$P{i}R'] = str(range_val)  # Range (OBLIGATOIRE)
            text_dict[f'$P{i}B'] = '32'  # Bits (OBLIGATOIRE)
            text_dict[f'$P{i}E'] = '0,0'  # Amplification type (OBLIGATOIRE)
            
            # Parametres OPTIONNELS mais RECOMMANDES
            text_dict[f'$P{i}G'] = '1'  # Gain
            text_dict[f'$P{i}D'] = 'Linear'  # Display type
            text_dict[f'$P{i}V'] = '1024'  # Voltage
            text_dict[f'$P{i}F'] = clean_name  # Filter name
            text_dict[f'$P{i}L'] = '0'  # Laser line
            text_dict[f'$P{i}O'] = 'PnDisplayName'  # Optical filter
            text_dict[f'$P{i}T'] = clean_name  # Detector type
        
        # Construire le TEXT segment avec encodage latin-1 (standard FCS)
        delimiter = '/'
        text_parts = []
        for key, value in text_dict.items():
            # Format: /key/value/
            text_parts.append(f'{key}{delimiter}{value}{delimiter}')
        
        text_str = delimiter + ''.join(text_parts)
        
        # Encoder en latin-1 (standard FCS)
        try:
            text_segment = text_str.encode('latin-1')
        except UnicodeEncodeError:
            # Si latin-1 echoue, remplacer les caracteres problematiques
            text_segment = text_str.encode('latin-1', errors='replace')
        
        # Padding pour aligner sur 8 bytes
        while len(text_segment) % 8 != 0:
            text_segment += b' '
        
        # Calculer les offsets
        header_size = 58
        text_start = header_size
        text_end = text_start + len(text_segment) - 1
        data_start = text_end + 1
        
        # Convertir les donnees en bytes (little endian float32)
        data_bytes = data.astype('<f4').tobytes()
        data_end = data_start + len(data_bytes) - 1
        
        # Ecrire le fichier FCS
        with open(filepath, 'wb') as f:
            # HEADER segment (58 bytes) - FCS 3.1
            header = f'FCS3.1    {text_start:8d}{text_end:8d}{data_start:8d}{data_end:8d}{0:8d}{0:8d}'
            # S'assurer que le header fait exactement 58 bytes
            header = header[:58].ljust(58)
            f.write(header.encode('ascii'))
            
            # TEXT segment
            f.write(text_segment)
            
            # DATA segment
            f.write(data_bytes)
        
        self._log(f"FCS ecrit (methode manuelle): {n_events} events, {n_channels} channels")
                
    def _log(self, message):
        if self.text_logs is not None:
            self.text_logs.append(message)
            scrollbar = self.text_logs.verticalScrollBar()
            scrollbar.setValue(scrollbar.maximum())
        # Toujours afficher dans la console pour debug
        print(f"[LOG] {message}")
        
    # =========================================================================
    # METHODES UMAP / t-SNE (avec Worker Thread)
    # =========================================================================
    
    def _start_embedding(self, method: str):
        """Lance le calcul d'embedding dans un thread separe."""
        if not self.result:
            self._fallback_visualization("Lancez d'abord une analyse FlowSOM")
            return
            
        # Verifier les dependances
        if method == 'umap' and not (UMAP_AVAILABLE or SCANPY_AVAILABLE):
            self._fallback_visualization("Installez umap-learn: pip install umap-learn")
            return
        if method == 'tsne' and not (TSNE_AVAILABLE or SCANPY_AVAILABLE):
            self._fallback_visualization("sklearn ou scanpy requis pour t-SNE")
            return
        
        # Arreter un calcul en cours
        if self.embedding_worker and self.embedding_worker.isRunning():
            self.embedding_worker.cancel()
            self.embedding_worker.terminate()
            self.embedding_worker.wait()
        
        # Recuperer le nombre de cellules configure
        if self.chk_all_cells.isChecked():
            # Utiliser toutes les cellules disponibles
            cell_data = self.result['cell_data']
            max_cells = cell_data.shape[0]
            self._log(f"ATTENTION: Utilisation de TOUTES les {max_cells:,} cellules (peut etre tres lent)")
        else:
            max_cells = self.spin_embedding_cells.value()
        
        # Afficher un message d'attente
        self.canvas.clear_figure()
        ax = self.canvas.fig.add_subplot(111)
        msg = f"Calcul {method.upper()} en cours ({max_cells:,} cellules)...\n"
        if max_cells > 20000:
            msg += "Attention: beaucoup de cellules, calcul tres long!"
        else:
            msg += "Ceci peut prendre 30-60 secondes."
        ax.text(0.5, 0.5, msg, 
               transform=ax.transAxes, ha='center', va='center', 
               color='#89b4fa', fontsize=14, fontweight='bold')
        ax.set_facecolor('#1e1e2e')
        self.canvas.fig.patch.set_facecolor('#1e1e2e')
        self.canvas.draw()
        
        # Lancer le worker
        cell_data = self.result['cell_data']
        cols_to_use = self.result['cols_to_use']
        
        self._log(f"Lancement {method.upper()} avec {max_cells:,} cellules...")
        
        self.embedding_worker = EmbeddingWorker(
            cell_data=cell_data,
            cols_to_use=cols_to_use,
            method=method,
            max_cells=max_cells,
            n_neighbors=15,
            min_dist=0.1,
            perplexity=30
        )
        self.embedding_worker.progress.connect(self._on_embedding_progress)
        self.embedding_worker.status.connect(self._on_embedding_status)
        self.embedding_worker.finished.connect(self._on_embedding_finished)
        self.embedding_worker.error.connect(self._on_embedding_error)
        self.embedding_worker.start()
        
    def _on_embedding_progress(self, value):
        """Callback pour la progression de l'embedding."""
        self.progress_bar.setValue(value)
        
    def _on_embedding_status(self, message):
        """Callback pour le statut de l'embedding."""
        self.lbl_status.setText(message)
        self._log(message)
        
    def _on_embedding_finished(self, result):
        """Callback quand l'embedding est termine."""
        self.embedding_result = result
        self.progress_bar.setValue(0)
        self._draw_embedding(result)
        
        # Stocker les coordonnees dans cell_data pour export
        if self.result and self.result.get('cell_data') is not None:
            cell_data = self.result['cell_data']
            method = result['method']
            embedding = result['embedding']
            indices = result['indices']
            
            # Creer un array de taille n_cells x 2 avec NaN par defaut
            n_cells = cell_data.shape[0]
            full_embedding = np.full((n_cells, 2), np.nan, dtype=np.float32)
            full_embedding[indices, :] = embedding
            
            # Stocker dans obsm
            key = f'X_{method}'
            cell_data.obsm[key] = full_embedding
            self._log(f"{method.upper()} coordonnees stockees pour export ({len(indices)} cellules)")
        
        self._log(f"{result['method'].upper()} termine avec succes!")
        
    def _on_embedding_error(self, error_message):
        """Callback en cas d'erreur d'embedding."""
        self._log(f"ERREUR embedding: {error_message}")
        self._fallback_visualization(error_message.split('\n')[0])
        self.progress_bar.setValue(0)
        
    def _draw_embedding(self, result):
        """Dessine le resultat de l'embedding."""
        self.canvas.clear_figure()
        
        embedding = result['embedding']
        metaclusters = result['metaclusters']
        method = result['method']
        
        ax = self.canvas.fig.add_subplot(111)
        
        # Palette de couleurs pour les clusters
        n_clusters = len(np.unique(metaclusters))
        colors = plt.cm.tab20(np.linspace(0, 1, max(n_clusters, 20)))
        
        # Scatter plot par metacluster
        for i in np.unique(metaclusters):
            mask = metaclusters == i
            ax.scatter(embedding[mask, 0], embedding[mask, 1], 
                      c=[colors[int(i) % 20]], label=f'MC {int(i)}',
                      s=5, alpha=0.6, edgecolors='none')
        
        ax.set_xlabel(f'{method.upper()}1', color='#cdd6f4', fontsize=12, fontweight='bold')
        ax.set_ylabel(f'{method.upper()}2', color='#cdd6f4', fontsize=12, fontweight='bold')
        ax.set_title(f'{method.upper()} - Metaclusters ({len(embedding):,} cellules)', 
                    fontsize=14, color='#cdd6f4', pad=20, fontweight='bold')
        ax.tick_params(colors='#cdd6f4', labelsize=10)
        ax.set_facecolor('#1e1e2e')
        ax.spines['bottom'].set_color('#45475a')
        ax.spines['left'].set_color('#45475a')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Ajuster les limites des axes avec une marge de 5%
        x_margin = (embedding[:, 0].max() - embedding[:, 0].min()) * 0.05
        y_margin = (embedding[:, 1].max() - embedding[:, 1].min()) * 0.05
        ax.set_xlim(embedding[:, 0].min() - x_margin, embedding[:, 0].max() + x_margin)
        ax.set_ylim(embedding[:, 1].min() - y_margin, embedding[:, 1].max() + y_margin)
        
        # Legende optimisée - en dehors du graphique ou compacte
        if n_clusters <= 15:
            legend = ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), 
                              fontsize=9, ncol=1,
                              facecolor='#313244', labelcolor='#cdd6f4', 
                              edgecolor='#45475a', markerscale=2.5,
                              framealpha=0.95)
        else:
            # Si trop de clusters, legende compacte en bas
            legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),
                              fontsize=7, ncol=min(n_clusters, 10),
                              facecolor='#313244', labelcolor='#cdd6f4', 
                              edgecolor='#45475a', markerscale=2,
                              framealpha=0.95)
        
        self.canvas.fig.patch.set_facecolor('#1e1e2e')
        self.canvas.fig.tight_layout(pad=2.0)
        self.canvas.draw()
        
    def _plot_umap(self):
        """Lance le calcul UMAP."""
        self._start_embedding('umap')
            
    def _plot_tsne(self):
        """Lance le calcul t-SNE."""
        self._start_embedding('tsne')
            
    def _embed_external_figure(self, fig, title):
        """Integre une figure matplotlib externe dans le canvas."""
        try:
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                tmp_path = tmp.name
            
            fig.savefig(tmp_path, dpi=150, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            plt.close(fig)
            
            img = plt.imread(tmp_path)
            ax = self.canvas.fig.add_subplot(111)
            ax.imshow(img)
            ax.axis('off')
            ax.set_title(title, fontsize=14, color='#cdd6f4', pad=15, fontweight='bold')
            self.canvas.fig.patch.set_facecolor('#1e1e2e')
            
            try:
                os.unlink(tmp_path)
            except:
                pass
        except Exception as e:
            self._fallback_visualization(str(e))
    
    # =========================================================================
    # METHODES SELECTION CLUSTERS INTERACTIFS
    # =========================================================================
    
    def _populate_cluster_list(self):
        """Remplit la liste des clusters apres l'analyse."""
        if not self.result or self.cluster_list is None:
            return
        
        # Bloquer les signaux pendant le remplissage pour eviter les appels en cascade
        self.cluster_list.blockSignals(True)
            
        n_clusters = self.result['n_clusters']
        self.cluster_list.clear()
        
        # Calculer les stats pour chaque cluster
        cell_data = self.result['cell_data']
        metaclustering = cell_data.obs['metaclustering'].values
        total_cells = len(metaclustering)
        
        for i in range(n_clusters):
            n_cells = (metaclustering == i).sum()
            pct = n_cells / total_cells * 100
            item = QListWidgetItem(f"Cluster {i} ({n_cells:,} - {pct:.1f}%)")
            item.setData(Qt.UserRole, i)
            self.cluster_list.addItem(item)
        
        # Reactiver les signaux
        self.cluster_list.blockSignals(False)
        
        # Ne pas selectionner automatiquement - laisser l'utilisateur choisir
        # Cela evite de bloquer l'interface avec les visualisations lourdes
        
    def _select_all_clusters(self):
        """Selectionne le premier cluster."""
        if self.cluster_list is None:
            return
        if self.cluster_list.count() > 0:
            self.cluster_list.setCurrentRow(0)
        
    def _deselect_all_clusters(self):
        """Deselectionne tous les clusters."""
        if self.cluster_list is None:
            return
        self.cluster_list.clearSelection()
        if self.lbl_cluster_info is not None:
            self.lbl_cluster_info.setText("Aucun cluster selectionne")
        self._clear_star_canvas()
        
    def _clear_star_canvas(self):
        """Efface le canvas du star plot."""
        self.star_canvas.clear_figure()
        ax = self.star_canvas.fig.add_subplot(111)
        ax.text(0.5, 0.5, "Selectionnez un cluster\npour afficher le Star Plot",
               transform=ax.transAxes, ha='center', va='center',
               color='#a6adc8', fontsize=11)
        ax.set_facecolor('#1e1e2e')
        self.star_canvas.fig.patch.set_facecolor('#1e1e2e')
        ax.axis('off')
        self.star_canvas.draw()
        
        # Effacer aussi le canvas FlowSOM Stars MST
        self.fsom_stars_canvas.clear_figure()
        ax2 = self.fsom_stars_canvas.fig.add_subplot(111)
        ax2.text(0.5, 0.5, "Le MST Star Plot FlowSOM apparaitra ici\napres selection d'un cluster",
               transform=ax2.transAxes, ha='center', va='center',
               color='#a6adc8', fontsize=11)
        ax2.set_facecolor('#1e1e2e')
        self.fsom_stars_canvas.fig.patch.set_facecolor('#1e1e2e')
        ax2.axis('off')
        self.fsom_stars_canvas.draw()
        
    def _on_cluster_selection_changed(self):
        """Callback quand la selection de clusters change."""
        # Activer le bouton MST si un cluster est selectionne
        if hasattr(self, 'btn_generate_mst') and self.btn_generate_mst is not None:
            has_selection = len(self.cluster_list.selectedItems()) > 0 if self.cluster_list else False
            self.btn_generate_mst.setEnabled(has_selection)
        
        # Appeler uniquement les visualisations legeres
        self._update_star_plot()
        self._update_cluster_comparison()
        # Ne PAS appeler _update_fsom_stars_mst automatiquement car c'est tres lourd
        # L'utilisateur peut le generer manuellement via le bouton
        self._show_mst_placeholder()
        
    def _update_star_plot(self):
        """Met a jour le Star Plot pour le cluster selectionne (version legere)."""
        if not self.result or self.cluster_list is None:
            return
            
        selected_items = self.cluster_list.selectedItems()
        if not selected_items:
            self._clear_star_canvas()
            return
            
        # Recuperer le cluster selectionne
        cluster_id = selected_items[0].data(Qt.UserRole)
        
        try:
            cell_data = self.result['cell_data']
            
            # Calculer les infos du cluster
            metaclustering = cell_data.obs['metaclustering'].values
            n_cells_cluster = (metaclustering == cluster_id).sum()
            total_cells = len(metaclustering)
            pct = n_cells_cluster / total_cells * 100
            
            self.lbl_cluster_info.setText(
                f"Cluster {cluster_id}: {n_cells_cluster:,} cellules ({pct:.2f}%)"
            )
            
            # Utiliser directement le graphique radar manuel (beaucoup plus leger)
            self.star_canvas.clear_figure()
            self._draw_manual_star_plot(cluster_id)
            
            self.star_canvas.fig.patch.set_facecolor('#1e1e2e')
            self.star_canvas.fig.tight_layout()
            self.star_canvas.draw()
            
        except Exception as e:
            self._log(f"Erreur Star Plot: {str(e)}")
            ax = self.star_canvas.fig.add_subplot(111)
            ax.text(0.5, 0.5, f"Erreur: {str(e)[:50]}...",
                   transform=ax.transAxes, ha='center', va='center',
                   color='#f38ba8', fontsize=10)
            ax.set_facecolor('#1e1e2e')
            self.star_canvas.fig.patch.set_facecolor('#1e1e2e')
            ax.axis('off')
            self.star_canvas.draw()
            
    def _draw_manual_star_plot(self, cluster_id: int):
        """Dessine un graphique radar (star plot) manuellement."""
        cell_data = self.result['cell_data']
        cols_to_use = self.result['cols_to_use']
        var_names = self.result.get('var_names', list(cell_data.var_names))
        
        # Recuperer les donnees du cluster
        metaclustering = cell_data.obs['metaclustering'].values
        mask = metaclustering == cluster_id
        
        X = cell_data.X
        if hasattr(X, 'toarray'):
            X = X.toarray()
        
        # Moyennes pour ce cluster
        cluster_means = np.nanmean(X[mask][:, cols_to_use], axis=0)
        
        # Normaliser
        cluster_means_norm = (cluster_means - np.nanmin(cluster_means)) / (np.nanmax(cluster_means) - np.nanmin(cluster_means) + 1e-10)
        
        # Noms des marqueurs
        markers = [var_names[i] for i in cols_to_use]
        
        # Limiter a 12 marqueurs max pour la lisibilite
        if len(markers) > 12:
            markers = markers[:12]
            cluster_means_norm = cluster_means_norm[:12]
        
        n_markers = len(markers)
        angles = np.linspace(0, 2 * np.pi, n_markers, endpoint=False).tolist()
        angles += angles[:1]  # Fermer le cercle
        
        values = cluster_means_norm.tolist()
        values += values[:1]  # Fermer le cercle
        
        ax = self.star_canvas.fig.add_subplot(111, polar=True)
        ax.fill(angles, values, color='#89b4fa', alpha=0.25)
        ax.plot(angles, values, 'o-', color='#89b4fa', linewidth=2, markersize=6)
        
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(markers, color='#cdd6f4', fontsize=8)
        ax.set_ylim(0, 1)
        ax.set_title(f'Expression - Cluster {cluster_id}', 
                    color='#cdd6f4', fontsize=12, fontweight='bold', pad=20)
        ax.tick_params(colors='#cdd6f4')
        ax.set_facecolor('#1e1e2e')
        ax.spines['polar'].set_color('#45475a')
        ax.grid(color='#45475a', alpha=0.3)
        
    def _update_cluster_comparison(self):
        """Met a jour le graphique de comparaison pour le cluster selectionne."""
        if not self.result or self.cluster_list is None or self.cluster_canvas is None:
            return
            
        selected_items = self.cluster_list.selectedItems()
        
        self.cluster_canvas.clear_figure()
        
        if not selected_items:
            ax = self.cluster_canvas.fig.add_subplot(111)
            ax.text(0.5, 0.5, "Selectionnez un cluster pour voir la comparaison",
                   transform=ax.transAxes, ha='center', va='center',
                   color='#a6adc8', fontsize=11)
            ax.set_facecolor('#1e1e2e')
            self.cluster_canvas.fig.patch.set_facecolor('#1e1e2e')
            self.cluster_canvas.draw()
            return
        
        cluster_id = selected_items[0].data(Qt.UserRole)
        
        try:
            cell_data = self.result['cell_data']
            cols_to_use = self.result['cols_to_use']
            var_names = self.result.get('var_names', list(cell_data.var_names))
            
            metaclustering = cell_data.obs['metaclustering'].values
            conditions = cell_data.obs['condition'].values
            
            mask_cluster = metaclustering == cluster_id
            
            X = cell_data.X
            if hasattr(X, 'toarray'):
                X = X.toarray()
            
            # Moyennes par condition pour ce cluster
            mask_healthy = (conditions == 'Sain') & mask_cluster
            mask_patho = (conditions == 'Pathologique') & mask_cluster
            
            markers = [var_names[i] for i in cols_to_use[:10]]  # Max 10 marqueurs
            
            healthy_means = np.nanmean(X[mask_healthy][:, cols_to_use[:10]], axis=0) if mask_healthy.sum() > 0 else np.zeros(len(cols_to_use[:10]))
            patho_means = np.nanmean(X[mask_patho][:, cols_to_use[:10]], axis=0) if mask_patho.sum() > 0 else np.zeros(len(cols_to_use[:10]))
            
            ax = self.cluster_canvas.fig.add_subplot(111)
            x = np.arange(len(markers))
            width = 0.35
            
            bars1 = ax.bar(x - width/2, healthy_means, width, label='NBM', 
                          color='#a6e3a1', edgecolor='white', linewidth=0.5)
            bars2 = ax.bar(x + width/2, patho_means, width, label='Patient', 
                          color='#f38ba8', edgecolor='white', linewidth=0.5)
            
            ax.set_xlabel('Marqueur', color='#cdd6f4', fontsize=10)
            ax.set_ylabel('Expression moyenne', color='#cdd6f4', fontsize=10)
            ax.set_title(f'Cluster {cluster_id} - Expression par Condition', 
                        fontsize=12, color='#cdd6f4', pad=10, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels(markers, color='#cdd6f4', fontsize=8, rotation=45, ha='right')
            ax.tick_params(colors='#cdd6f4')
            ax.legend(facecolor='#313244', labelcolor='#cdd6f4', edgecolor='#45475a', fontsize=8)
            ax.set_facecolor('#1e1e2e')
            ax.spines['bottom'].set_color('#45475a')
            ax.spines['left'].set_color('#45475a')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            self.cluster_canvas.fig.patch.set_facecolor('#1e1e2e')
            self.cluster_canvas.fig.tight_layout()
            
        except Exception as e:
            ax = self.cluster_canvas.fig.add_subplot(111)
            ax.text(0.5, 0.5, f"Erreur: {str(e)}", 
                   transform=ax.transAxes, ha='center', va='center',
                   color='#f38ba8', fontsize=10)
            ax.set_facecolor('#1e1e2e')
            self.cluster_canvas.fig.patch.set_facecolor('#1e1e2e')
            
        self.cluster_canvas.draw()
    
    def _show_mst_placeholder(self):
        """Affiche un placeholder pour le MST (generation a la demande)."""
        if self.fsom_stars_canvas is None:
            return
        
        self.fsom_stars_canvas.clear_figure()
        ax = self.fsom_stars_canvas.fig.add_subplot(111)
        
        selected_items = self.cluster_list.selectedItems() if self.cluster_list else []
        if selected_items:
            cluster_id = selected_items[0].data(Qt.UserRole)
            ax.text(0.5, 0.5, f"MST pour Cluster {cluster_id}\n\nCliquez sur 'Generer MST' pour afficher",
                   transform=ax.transAxes, ha='center', va='center',
                   color='#a6adc8', fontsize=11)
        else:
            ax.text(0.5, 0.5, "Selectionnez un cluster",
                   transform=ax.transAxes, ha='center', va='center',
                   color='#a6adc8', fontsize=11)
        
        ax.set_facecolor('#1e1e2e')
        self.fsom_stars_canvas.fig.patch.set_facecolor('#1e1e2e')
        ax.axis('off')
        self.fsom_stars_canvas.draw()

    def _update_fsom_stars_mst(self):
        """
        Met a jour l'image FlowSOM Stars MST pour le cluster selectionne.
        Utilise exactement:
            fsom_subset = fsom.subset(fsom.get_cell_data().obs["metaclustering"] == cluster_id)
            p = fs.pl.plot_stars(fsom_subset, background_values=fsom_subset.get_cluster_data().obs.metaclustering)
        """
        if not self.result or self.cluster_list is None or self.fsom_stars_canvas is None:
            return
            
        selected_items = self.cluster_list.selectedItems()
        
        self.fsom_stars_canvas.clear_figure()
        
        if not selected_items:
            ax = self.fsom_stars_canvas.fig.add_subplot(111)
            ax.text(0.5, 0.5, "Selectionnez un cluster pour voir le MST",
                   transform=ax.transAxes, ha='center', va='center',
                   color='#a6adc8', fontsize=11)
            ax.set_facecolor('#1e1e2e')
            self.fsom_stars_canvas.fig.patch.set_facecolor('#1e1e2e')
            ax.axis('off')
            self.fsom_stars_canvas.draw()
            return
        
        cluster_id = selected_items[0].data(Qt.UserRole)
        
        try:
            fsom = self.result['fsom']
            
            # Exactement comme dans l'exemple:
            # fsom_subset = fsom.subset(fsom.get_cell_data().obs["metaclustering"] == cluster_id)
            fsom_subset = fsom.subset(fsom.get_cell_data().obs["metaclustering"] == cluster_id)
            
            # p = fs.pl.plot_stars(fsom_subset, background_values=fsom_subset.get_cluster_data().obs.metaclustering)
            fig = fs.pl.plot_stars(
                fsom_subset, 
                background_values=fsom_subset.get_cluster_data().obs.metaclustering
            )
            
            # Sauvegarder temporairement l'image PNG
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                tmp_path = tmp.name
                
                # Si fig est une figure matplotlib
                if hasattr(fig, 'savefig'):
                    fig.savefig(tmp_path, dpi=150, bbox_inches='tight', 
                               facecolor='white', edgecolor='none')
                    plt.close(fig)
                else:
                    # Si c'est retourne autrement, essayer plt.gcf()
                    plt.savefig(tmp_path, dpi=150, bbox_inches='tight',
                               facecolor='white', edgecolor='none')
                    plt.close()
            
            # Charger l'image dans le canvas
            img = plt.imread(tmp_path)
            ax = self.fsom_stars_canvas.fig.add_subplot(111)
            ax.imshow(img)
            ax.axis('off')
            ax.set_title(f'FlowSOM MST - Cluster {cluster_id}', 
                        color='#cdd6f4', fontsize=14, fontweight='bold', pad=10)
            
            self.fsom_stars_canvas.fig.patch.set_facecolor('#1e1e2e')
            self.fsom_stars_canvas.fig.tight_layout()
            
            # Nettoyer le fichier temporaire
            try:
                os.unlink(tmp_path)
            except:
                pass
                
            self._log(f"MST Star Plot genere pour Cluster {cluster_id}")
            
        except Exception as e:
            self._log(f"Erreur MST Star Plot: {str(e)}")
            ax = self.fsom_stars_canvas.fig.add_subplot(111)
            ax.text(0.5, 0.5, f"Erreur generation MST:\n{str(e)[:60]}...",
                   transform=ax.transAxes, ha='center', va='center',
                   color='#f38ba8', fontsize=10, wrap=True)
            ax.set_facecolor('#1e1e2e')
            self.fsom_stars_canvas.fig.patch.set_facecolor('#1e1e2e')
            ax.axis('off')
            
        self.fsom_stars_canvas.draw()
    
    # =========================================================================
    # METHODES TRANSFORMATION DES DONNEES
    # =========================================================================
    
    def _apply_transformation(self, X: np.ndarray) -> np.ndarray:
        """Applique la transformation selectionnee aux donnees."""
        transform_type = self.combo_transform.currentText()
        
        if transform_type == "Aucune":
            return X
        elif "Arcsinh" in transform_type:
            if "cofactor=5" in transform_type:
                cofactor = 5.0
            elif "cofactor=150" in transform_type:
                cofactor = 150.0
            else:
                cofactor = self.spin_cofactor.value()
            self._log(f"Application Arcsinh (cofactor={cofactor})")
            return DataTransformer.arcsinh_transform(X, cofactor)
        elif transform_type == "Logicle":
            self._log("Application transformation Logicle")
            return DataTransformer.logicle_transform(X)
        elif transform_type == "Log10":
            self._log("Application transformation Log10")
            return DataTransformer.log_transform(X)
        elif transform_type == "Z-score":
            self._log("Application normalisation Z-score")
            return DataTransformer.zscore_normalize(X)
        else:
            return X
    
    # =========================================================================
    # METHODES SUIVI PATIENT
    # =========================================================================
    
    def _save_patient_timepoint(self):
        """Sauvegarde le timepoint actuel pour le patient."""
        patient_id = self.edit_patient_id.text().strip()
        timepoint = self.edit_timepoint.text().strip()
        
        if not patient_id:
            QMessageBox.warning(self, "Attention", "Entrez un ID patient.")
            return
        if not timepoint:
            QMessageBox.warning(self, "Attention", "Entrez un nom de timepoint.")
            return
        if not self.result:
            QMessageBox.warning(self, "Attention", "Lancez d'abord une analyse.")
            return
        
        # Creer ou recuperer le patient
        if patient_id not in self.patients_db:
            self.patients_db[patient_id] = PatientData(patient_id)
        
        patient = self.patients_db[patient_id]
        self.current_patient = patient
        
        # Ajouter le timepoint
        files = self.healthy_files + self.pathological_files
        patient.add_timepoint(timepoint, files, datetime.now().strftime("%Y-%m-%d"))
        
        # Enregistrer les resultats
        results_summary = {
            'n_clusters': self.result.get('n_clusters'),
            'n_cells': self.result['cell_data'].shape[0] if self.result.get('cell_data') is not None else 0,
            'lsc_percentage': self.lsc_result.get('lsc_percentage', 0) if self.lsc_result else 0
        }
        patient.set_timepoint_results(timepoint, results_summary)
        patient.timepoints[timepoint]['lsc_score'] = self.lsc_result.get('lsc_percentage', 0) if self.lsc_result else None
        
        # Sauvegarder sur disque
        self._save_patients_to_disk()
        
        self._log(f"Timepoint {timepoint} sauvegarde pour patient {patient_id}")
        # Mettre a jour la timeline si current_patient existe
        if self.current_patient:
            self._update_patient_timeline({'timepoints': list(self.current_patient.timepoints.values()), 'patient_id': patient_id})
        
        QMessageBox.information(self, "Succes", f"Timepoint {timepoint} sauvegarde!")
        
    def _load_patient_history(self):
        """Charge l'historique d'un patient depuis un fichier."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Charger historique patient", "",
            "Fichiers JSON (*.json);;Tous les fichiers (*)"
        )
        
        if file_path:
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                
                patient = PatientData.from_dict(data)
                self.patients_db[patient.patient_id] = patient
                self.current_patient = patient
                self.edit_patient_id.setText(patient.patient_id)
                
                # Convertir les timepoints pour la mise a jour
                tp_list = [{'id': tid, 'date': tp.get('date', ''), 'files_analyzed': tp.get('files', []), 'n_clusters': tp.get('results', {}).get('n_clusters', 0), 'notes': tp.get('notes', '')} for tid, tp in patient.timepoints.items()]
                self._update_patient_timeline({'timepoints': tp_list, 'patient_id': patient.patient_id})
                self._log(f"Historique patient {patient.patient_id} charge ({len(patient.timepoints)} timepoints)")
                
            except Exception as e:
                QMessageBox.critical(self, "Erreur", f"Erreur chargement: {str(e)}")
                
    def _save_patients_to_disk(self):
        """Sauvegarde tous les patients sur disque."""
        data_folder = self._get_patient_data_folder()
        for patient_id, patient in self.patients_db.items():
            try:
                filename = f"{patient_id.replace(' ', '_')}.json"
                filepath = data_folder / filename
                
                with open(filepath, 'w', encoding='utf-8') as f:
                    json.dump(patient.to_dict(), f, indent=2, ensure_ascii=False)
            except Exception as e:
                self._log(f"Erreur sauvegarde patient {patient_id}: {str(e)}")
                
    def _update_patient_timeline_from_current(self):
        """Met a jour l'affichage de la timeline patient depuis current_patient."""
        if not self.current_patient or self.patient_timeline is None:
            return
        
        patient = self.current_patient
        timeline = patient.get_timeline()
        
        self.patient_timeline.setRowCount(len(timeline))
        
        lsc_values = []
        dates = []
        
        for i, (tid, date) in enumerate(timeline):
            tp = patient.timepoints[tid]
            lsc = tp.get('lsc_score')
            
            self.patient_timeline.setItem(i, 0, QTableWidgetItem(tid))
            self.patient_timeline.setItem(i, 1, QTableWidgetItem(date))
            self.patient_timeline.setItem(i, 2, QTableWidgetItem(f"{lsc:.2f}%" if lsc else "--"))
            self.patient_timeline.setItem(i, 3, QTableWidgetItem(tp.get('notes', '')))
            
            if lsc is not None:
                lsc_values.append(lsc)
                dates.append(tid)
        
        # Dessiner l'evolution temporelle
        self._draw_patient_evolution(dates, lsc_values)
        
    def _draw_patient_evolution(self, dates: List[str], lsc_values: List[float]):
        """Dessine l'evolution du score LSC dans le temps."""
        if self.patient_canvas is None:
            return
        self.patient_canvas.clear_figure()
        
        if not dates or not lsc_values:
            ax = self.patient_canvas.fig.add_subplot(111)
            ax.text(0.5, 0.5, "Aucune donnee temporelle disponible",
                   transform=ax.transAxes, ha='center', va='center',
                   color='#a6adc8', fontsize=11)
            ax.set_facecolor('#1e1e2e')
            self.patient_canvas.fig.patch.set_facecolor('#1e1e2e')
            self.patient_canvas.draw()
            return
        
        ax = self.patient_canvas.fig.add_subplot(111)
        
        x = np.arange(len(dates))
        ax.plot(x, lsc_values, 'o-', color='#f38ba8', linewidth=2, markersize=10)
        ax.fill_between(x, lsc_values, alpha=0.2, color='#f38ba8')
        
        # Ajouter les valeurs
        for i, val in enumerate(lsc_values):
            ax.annotate(f'{val:.2f}%', (i, val), textcoords="offset points",
                       xytext=(0, 10), ha='center', color='#cdd6f4', fontsize=10)
        
        # Ligne de seuil
        ax.axhline(y=0.1, color='#a6e3a1', linestyle='--', alpha=0.7, label='Seuil MRD negatif')
        ax.axhline(y=1.0, color='#f9e2af', linestyle='--', alpha=0.7, label='Seuil alerte')
        
        ax.set_xlabel('Timepoint', color='#cdd6f4', fontsize=11)
        ax.set_ylabel('Score LSC (%)', color='#cdd6f4', fontsize=11)
        ax.set_title(f'Evolution LSC - Patient {self.current_patient.patient_id}',
                    color='#cdd6f4', fontsize=14, fontweight='bold', pad=15)
        ax.set_xticks(x)
        ax.set_xticklabels(dates, color='#cdd6f4')
        ax.tick_params(colors='#cdd6f4')
        ax.legend(facecolor='#313244', labelcolor='#cdd6f4', edgecolor='#45475a', loc='upper right')
        ax.set_facecolor('#1e1e2e')
        ax.spines['bottom'].set_color('#45475a')
        ax.spines['left'].set_color('#45475a')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.2, color='#45475a')
        
        self.patient_canvas.fig.patch.set_facecolor('#1e1e2e')
        self.patient_canvas.fig.tight_layout()
        self.patient_canvas.draw()
    
    # =========================================================================
    # METHODES EXPORT PDF
    # =========================================================================
    
    def _export_pdf_report(self):
        """Exporte un rapport PDF complet."""
        if not REPORTLAB_AVAILABLE:
            QMessageBox.warning(self, "Attention", 
                "ReportLab n'est pas installe.\npip install reportlab")
            return
        
        if not self.result:
            QMessageBox.warning(self, "Attention", "Lancez d'abord une analyse.")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Exporter Rapport PDF", 
            f"FlowSOM_Report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf",
            "Fichiers PDF (*.pdf)"
        )
        
        if not file_path:
            return
        
        try:
            self._log("Generation du rapport PDF...")
            
            # Sauvegarder les figures temporairement
            fig_paths = []
            
            # Figure principale
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                self.canvas.fig.savefig(tmp.name, dpi=800, bbox_inches='tight',
                                       facecolor='white', edgecolor='none')
                fig_paths.append(tmp.name)
            
            # Figure LSC si disponible
            if self.lsc_result and hasattr(self, 'lsc_canvas') and self.lsc_canvas:
                with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                    self.lsc_canvas.fig.savefig(tmp.name, dpi=800, bbox_inches='tight',
                                               facecolor='white', edgecolor='none')
                    fig_paths.append(tmp.name)
            
            # Infos patient
            patient_info = None
            if self.edit_patient_id and self.edit_patient_id.text().strip():
                patient_info = {
                    'patient_id': self.edit_patient_id.text().strip(),
                    'date': datetime.now().strftime('%d/%m/%Y'),
                    'notes': self.edit_timepoint.text().strip() if self.edit_timepoint else ''
                }
            
            # Generer le PDF
            generator = PDFReportGenerator(file_path)
            success = generator.generate_report(
                self.result,
                fig_paths,
                patient_info,
                self.lsc_result
            )
            
            # Nettoyer les fichiers temporaires
            for path in fig_paths:
                try:
                    os.unlink(path)
                except:
                    pass
            
            if success:
                self._log(f"Rapport PDF genere: {file_path}")
                QMessageBox.information(self, "Succes", f"Rapport PDF exporte:\n{file_path}")
            else:
                QMessageBox.warning(self, "Erreur", "Erreur lors de la generation du PDF.")
                
        except Exception as e:
            self._log(f"Erreur export PDF: {str(e)}")
            QMessageBox.critical(self, "Erreur", f"Erreur: {str(e)}")
    
    # =========================================================================
    # METHODES ELN CONFORMITE - NBM REFERENCE [2]
    # =========================================================================
    
    def _load_nbm_folder(self):
        """Charge un dossier de fichiers NBM pour construire la reference."""
        folder = QFileDialog.getExistingDirectory(
            self, "Selectionner le dossier NBM (Normal Bone Marrow)", "",
            QFileDialog.ShowDirsOnly
        )
        if folder:
            self.nbm_files = self._get_fcs_files(folder)
            n = len(self.nbm_files)
            
            min_required = self.config['quality_control'].get('nbm_samples_required', 15)
            
            if self.lbl_nbm_status:
                if n >= min_required:
                    self.lbl_nbm_status.setText(f"NBM: {n} fichiers (OK)")
                    self.lbl_nbm_status.setStyleSheet("color: #a6e3a1; font-size: 9pt; padding: 5px;")
                else:
                    self.lbl_nbm_status.setText(f"NBM: {n}/{min_required} fichiers (insuffisant)")
                    self.lbl_nbm_status.setStyleSheet("color: #f9e2af; font-size: 9pt; padding: 5px;")
            
            if self.btn_build_nbm_ref:
                self.btn_build_nbm_ref.setEnabled(n > 0)
            
            self._log(f"Dossier NBM: {folder} ({n} fichiers)")
    
    def _build_nbm_reference(self):
        """Construit la reference NBM frozen MST."""
        if not self.nbm_files:
            QMessageBox.warning(self, "Attention", "Chargez d'abord un dossier NBM.")
            return
        
        if not FLOWSOM_AVAILABLE:
            QMessageBox.critical(self, "Erreur", "FlowSOM n'est pas installe.")
            return
        
        try:
            self._log("Construction de la reference NBM...")
            self.statusBar().showMessage("Construction reference NBM...")
            QApplication.processEvents()
            
            # Charger les donnees NBM
            nbm_data_list = []
            var_names = None
            
            for fpath in self.nbm_files:
                try:
                    adata = fs.io.read_FCS(fpath)
                    nbm_data_list.append(adata.X)
                    if var_names is None:
                        var_names = list(adata.var_names)
                except Exception as e:
                    self._log(f"Erreur lecture {Path(fpath).name}: {e}")
            
            if not nbm_data_list:
                QMessageBox.critical(self, "Erreur", "Aucun fichier NBM valide.")
                return
            
            # Construire la reference
            self.nbm_reference = NBMReferenceBuilder(n_nodes=100, seed=42)
            success = self.nbm_reference.build_reference(nbm_data_list, var_names)
            
            if success:
                if self.lbl_nbm_status:
                    self.lbl_nbm_status.setText(f"Reference NBM: Construite ({len(nbm_data_list)} samples)")
                    self.lbl_nbm_status.setStyleSheet("color: #a6e3a1; font-weight: bold; font-size: 9pt; padding: 5px;")
                
                # Initialiser le clone tracker avec la reference
                self.clone_tracker = ClonePersistenceTracker(
                    nbm_reference=self.nbm_reference,
                    mrd_detector=self.mrd_detector
                )
                
                self._log(f"Reference NBM construite avec {len(nbm_data_list)} echantillons")
                self.statusBar().showMessage("Reference NBM construite avec succes!")
                
                # Proposer de sauvegarder
                reply = QMessageBox.question(
                    self, "Sauvegarder Reference",
                    "Voulez-vous sauvegarder la reference NBM pour reutilisation?",
                    QMessageBox.Yes | QMessageBox.No
                )
                
                if reply == QMessageBox.Yes:
                    file_path, _ = QFileDialog.getSaveFileName(
                        self, "Sauvegarder Reference NBM",
                        "nbm_reference.json",
                        "Fichiers JSON (*.json)"
                    )
                    if file_path:
                        self.nbm_reference.save_reference(file_path)
                        self._log(f"Reference NBM sauvegardee: {file_path}")
            else:
                QMessageBox.critical(self, "Erreur", "Erreur construction reference NBM.")
                
        except Exception as e:
            self._log(f"Erreur construction NBM: {str(e)}")
            QMessageBox.critical(self, "Erreur", f"Erreur: {str(e)}")
    
    def _load_nbm_reference(self):
        """Charge une reference NBM existante."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Charger Reference NBM",
            "",
            "Fichiers JSON (*.json)"
        )
        
        if not file_path:
            return
        
        try:
            self.nbm_reference = NBMReferenceBuilder()
            success = self.nbm_reference.load_reference(file_path)
            
            if success:
                if self.lbl_nbm_status:
                    self.lbl_nbm_status.setText(
                        f"Reference NBM: Chargee ({self.nbm_reference.nbm_samples_count} samples)")
                    self.lbl_nbm_status.setStyleSheet("color: #a6e3a1; font-weight: bold; font-size: 9pt; padding: 5px;")
                
                # Initialiser le clone tracker
                self.clone_tracker = ClonePersistenceTracker(
                    nbm_reference=self.nbm_reference,
                    mrd_detector=self.mrd_detector
                )
                
                self._log(f"Reference NBM chargee depuis {file_path}")
            else:
                QMessageBox.critical(self, "Erreur", "Erreur chargement reference NBM.")
                
        except Exception as e:
            self._log(f"Erreur chargement NBM: {str(e)}")
            QMessageBox.critical(self, "Erreur", f"Erreur: {str(e)}")
    
    # =========================================================================
    # METHODES ELN CONFORMITE - CLONE PERSISTENCE [2]
    # =========================================================================
    
    def _analyze_clone_persistence(self):
        """Analyse la persistance clonale entre T0 et FU."""
        if not self.result:
            QMessageBox.warning(self, "Attention", "Lancez d'abord une analyse.")
            return
        
        if not self.clone_tracker:
            QMessageBox.warning(self, "Attention", 
                "Construisez d'abord une reference NBM pour l'analyse de persistance clonale.")
            return
        
        if not self.edit_patient_id or not self.edit_timepoint:
            return
        
        patient_id = self.edit_patient_id.text().strip()
        timepoint_id = self.edit_timepoint.text().strip()
        
        if not patient_id or not timepoint_id:
            QMessageBox.warning(self, "Attention", "Renseignez l'ID patient et le timepoint.")
            return
        
        try:
            self._log(f"Analyse persistance clonale pour {patient_id}...")
            
            # Obtenir les frequences par node pour ce timepoint
            fsom = self.result.get('fsom')
            if fsom is None:
                return
            
            cell_data = fsom.get_cell_data()
            node_assignments = cell_data.obs['clustering'].values
            total_cells = len(node_assignments)
            
            node_frequencies = {}
            for node_id in range(100):  # 100 nodes standard
                count = np.sum(node_assignments == str(node_id))
                node_frequencies[node_id] = count / total_cells
            
            # Ajouter au tracker
            self.clone_tracker.add_timepoint(timepoint_id, node_frequencies)
            
            # Si on a T0 et un FU, faire l'analyse
            timepoints = list(self.clone_tracker.timepoints.keys())
            if len(timepoints) >= 2:
                # Comparer avec T0 (premier timepoint)
                t0_id = sorted(timepoints)[0]
                fu_id = timepoint_id
                
                if t0_id != fu_id:
                    result = self.clone_tracker.track_clone_persistence(t0_id, fu_id)
                    
                    # Afficher les resultats
                    self._display_clone_persistence_results(result)
                    
                    # Mettre a jour le label MRD
                    if self.lbl_mrd_result:
                        mrd_status = result.get('mrd_status', 'Unknown')
                        mrd_percent = result.get('total_mrd_percent', 0)
                        
                        if mrd_status == 'Undetectable':
                            color = "#a6e3a1"
                        elif mrd_status == 'Below LOQ':
                            color = "#f9e2af"
                        else:
                            color = "#f38ba8"
                        
                        self.lbl_mrd_result.setText(f"MRD: {mrd_status} ({mrd_percent:.4f}%)")
                        self.lbl_mrd_result.setStyleSheet(f"""
                            font-size: 11pt;
                            font-weight: bold;
                            color: {color};
                            background: rgba({color[1:]}, 0.1);
                            padding: 8px;
                            border-radius: 6px;
                        """)
            else:
                self._log(f"Timepoint {timepoint_id} enregistre. Ajoutez un 2eme timepoint pour l'analyse.")
                QMessageBox.information(self, "Info", 
                    f"Timepoint '{timepoint_id}' enregistre.\n"
                    "Ajoutez un deuxieme timepoint pour l'analyse de persistance.")
                
        except Exception as e:
            self._log(f"Erreur analyse persistance: {str(e)}")
            QMessageBox.critical(self, "Erreur", f"Erreur: {str(e)}")
    
    def _display_clone_persistence_results(self, result: Dict[str, Any]):
        """Affiche les resultats de persistance clonale."""
        self._log("=" * 50)
        self._log("ANALYSE PERSISTANCE CLONALE")
        self._log("=" * 50)
        self._log(f"Diagnostic: {result['diagnosis_id']} -> Follow-up: {result['followup_id']}")
        self._log(f"MRD Total: {result['total_mrd_percent']:.4f}%")
        self._log(f"Status: {result['mrd_status']}")
        self._log("")
        
        # Clones persistants
        self._log(f"Clones persistants: {len(result['persistent_nodes'])}")
        for node in result['persistent_nodes'][:5]:  # Top 5
            if node['mrd_positive']:
                self._log(f"  Node {node['node_id']}: {node['followup_freq']:.4%} "
                         f"(fold-change NBM: {node['fold_change_to_nbm']:.1f}x) - MRD+")
        
        # Clones resolus
        self._log(f"Clones resolus: {len(result['resolved_nodes'])}")
        
        # Nouveaux clones
        if result['new_nodes']:
            self._log(f"ATTENTION - Nouveaux clones detectes: {len(result['new_nodes'])}")
            for node in result['new_nodes'][:3]:
                self._log(f"  Node {node['node_id']}: {node['followup_freq']:.4%}")
        
        self._log("=" * 50)


# =============================================================================
# POINT D'ENTREE
# =============================================================================

def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    
    # Palette sombre
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(30, 30, 46))
    palette.setColor(QPalette.WindowText, QColor(205, 214, 244))
    palette.setColor(QPalette.Base, QColor(49, 50, 68))
    palette.setColor(QPalette.AlternateBase, QColor(69, 71, 90))
    palette.setColor(QPalette.ToolTipBase, QColor(49, 50, 68))
    palette.setColor(QPalette.ToolTipText, QColor(205, 214, 244))
    palette.setColor(QPalette.Text, QColor(205, 214, 244))
    palette.setColor(QPalette.Button, QColor(69, 71, 90))
    palette.setColor(QPalette.ButtonText, QColor(205, 214, 244))
    palette.setColor(QPalette.Highlight, QColor(137, 180, 250))
    palette.setColor(QPalette.HighlightedText, QColor(30, 30, 46))
    app.setPalette(palette)
    
    window = FlowSOMApp()
    window.show()
    
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
