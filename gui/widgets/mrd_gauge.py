# -*- coding: utf-8 -*-
"""
mrd_gauge.py — Widget gauge MRD par méthode.

Affiche pour une méthode donnée :
  - Nom de la méthode
  - Pourcentage MRD (barre de progression + valeur numérique)
  - Badge POSITIF / NÉGATIF / INDÉTECTABLE
  - Nombre de nœuds MRD et cellules impliquées
"""
from __future__ import annotations

from typing import Any, Dict, Optional

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QProgressBar, QFrame,
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QColor


# Couleurs Catppuccin Mocha (cohérence avec styles.py)
_RED = "#f38ba8"
_YELLOW = "#f9e2af"
_GREEN = "#a6e3a1"
_BLUE = "#89b4fa"
_SUBTEXT = "#a6adc8"
_SURFACE0 = "#313244"
_SURFACE1 = "#45475a"
_BASE = "#1e1e2e"
_MANTLE = "#181825"
_TEXT = "#cdd6f4"


class MRDGauge(QWidget):
    """
    Widget compact affichant le résultat MRD pour une méthode.

    Usage :
        gauge = MRDGauge(method_name="ELN 2025")
        gauge.update_data({
            "pct": 0.42,
            "n_cells": 1234,
            "n_nodes": 3,
            "positive": True,
            "positivity_threshold": 0.1,  # optionnel
        })
    """

    def __init__(self, method_name: str, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.method_name = method_name
        self._build_ui()
        self._set_empty()

    def _build_ui(self) -> None:
        self.setObjectName("mrdGaugeCard")
        self._apply_card_style(status="waiting")

        root = QVBoxLayout(self)
        root.setContentsMargins(18, 16, 18, 16)
        root.setSpacing(10)

        # ── Titre méthode ──
        self.lbl_method = QLabel(self.method_name)
        self.lbl_method.setFont(QFont("Segoe UI", 10, QFont.Bold))
        self.lbl_method.setStyleSheet(
            f"color: #7a9ccc; background: transparent; letter-spacing: 0.12em; text-transform: uppercase;"
        )
        self.lbl_method.setAlignment(Qt.AlignCenter)
        root.addWidget(self.lbl_method)

        # ── Badge statut ──
        self.lbl_badge = QLabel("EN ATTENTE")
        self.lbl_badge.setAlignment(Qt.AlignCenter)
        self.lbl_badge.setFont(QFont("Segoe UI", 8, QFont.Bold))
        self.lbl_badge.setFixedHeight(22)
        self.lbl_badge.setStyleSheet(
            f"color: {_SUBTEXT}; background: rgba(69,71,90,0.6); "
            f"border-radius: 5px; padding: 0 10px; letter-spacing: 0.08em;"
        )
        root.addWidget(self.lbl_badge)

        # ── Valeur MRD ──
        self.lbl_pct = QLabel("—")
        self.lbl_pct.setAlignment(Qt.AlignCenter)
        self.lbl_pct.setFont(QFont("Segoe UI", 32, QFont.Bold))
        self.lbl_pct.setStyleSheet(f"color: {_TEXT}; background: transparent; letter-spacing: -0.02em;")
        root.addWidget(self.lbl_pct)

        # ── Barre de progression ──
        self.progress = QProgressBar()
        self.progress.setRange(0, 1000)
        self.progress.setValue(0)
        self.progress.setTextVisible(False)
        self.progress.setFixedHeight(10)
        self.progress.setStyleSheet(f"""
            QProgressBar {{
                background: rgba(30, 32, 50, 0.9);
                border-radius: 5px;
                border: none;
            }}
            QProgressBar::chunk {{
                border-radius: 5px;
                background: {_BLUE};
            }}
        """)
        root.addWidget(self.progress)

        # ── Infos secondaires ──
        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setStyleSheet("color: rgba(137,180,250,0.08); max-height: 1px;")
        root.addWidget(sep)

        info_row = QHBoxLayout()
        info_row.setSpacing(12)

        self.lbl_nodes = QLabel("— nœuds")
        self.lbl_nodes.setStyleSheet(
            f"color: #4a4c70; background: transparent; font-size: 10px; font-weight: 600;"
        )
        self.lbl_nodes.setAlignment(Qt.AlignCenter)

        self.lbl_cells = QLabel("— cellules")
        self.lbl_cells.setStyleSheet(
            f"color: #4a4c70; background: transparent; font-size: 10px; font-weight: 600;"
        )
        self.lbl_cells.setAlignment(Qt.AlignCenter)

        info_row.addWidget(self.lbl_nodes)
        info_row.addWidget(self.lbl_cells)
        root.addLayout(info_row)

    def _apply_card_style(self, status: str = "waiting") -> None:
        if status == "positive":
            border_color = "rgba(243, 139, 168, 0.5)"
            bg_top = "rgba(50, 28, 36, 0.85)"
            bg_bot = "rgba(38, 20, 28, 0.85)"
            top_border = "rgba(243, 139, 168, 0.7)"
        elif status == "negative":
            border_color = "rgba(166, 227, 161, 0.4)"
            bg_top = "rgba(28, 50, 30, 0.8)"
            bg_bot = "rgba(20, 38, 22, 0.8)"
            top_border = "rgba(166, 227, 161, 0.55)"
        elif status == "low":
            border_color = "rgba(249, 226, 175, 0.4)"
            bg_top = "rgba(50, 44, 20, 0.8)"
            bg_bot = "rgba(38, 33, 14, 0.8)"
            top_border = "rgba(249, 226, 175, 0.55)"
        else:  # waiting
            border_color = "rgba(69, 71, 90, 0.6)"
            bg_top = "rgba(44, 46, 64, 0.7)"
            bg_bot = "rgba(32, 34, 52, 0.7)"
            top_border = "rgba(80, 82, 102, 0.5)"
        self.setStyleSheet(f"""
            QWidget#mrdGaugeCard {{
                background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                    stop:0 {bg_top}, stop:1 {bg_bot});
                border-radius: 12px;
                border: 1px solid {border_color};
                border-top: 2px solid {top_border};
            }}
        """)

    # ------------------------------------------------------------------

    def update_data(self, data: Dict[str, Any]) -> None:
        """Met à jour l'affichage avec les données MRD de la méthode."""
        pct = data.get("pct", 0.0)
        n_nodes = data.get("n_nodes", 0)
        n_cells = data.get("n_cells", 0)
        positive = data.get("positive", False)
        low_level = data.get("low_level", False)
        threshold = data.get("positivity_threshold", None)

        # Valeur numérique
        self.lbl_pct.setText(f"{pct:.4f} %")

        # Barre (max représente 5%, clampé)
        bar_val = min(int(pct * 200), 1000)  # 5% → 1000
        self.progress.setValue(bar_val)

        # Badge et couleurs
        if positive or n_nodes > 0:
            badge_text = "POSITIF"
            badge_color = _RED
            badge_bg = "rgba(243, 139, 168, 0.22)"
            badge_border = "rgba(243, 139, 168, 0.5)"
            bar_color = f"qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 {_RED}, stop:1 #eb6f92)"
            pct_color = "#f7a8bf"
            card_status = "positive"
        elif low_level:
            badge_text = "INDÉTECTABLE"
            badge_color = _YELLOW
            badge_bg = "rgba(249, 226, 175, 0.18)"
            badge_border = "rgba(249, 226, 175, 0.4)"
            bar_color = f"qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 {_YELLOW}, stop:1 #f0c060)"
            pct_color = _YELLOW
            card_status = "low"
        else:
            badge_text = "NÉGATIF"
            badge_color = _GREEN
            badge_bg = "rgba(166, 227, 161, 0.18)"
            badge_border = "rgba(166, 227, 161, 0.4)"
            bar_color = f"qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 {_GREEN}, stop:1 #78c878)"
            pct_color = "#b8f0b4"
            card_status = "negative"

        self._apply_card_style(card_status)
        self.lbl_pct.setStyleSheet(f"color: {pct_color}; background: transparent; letter-spacing: -0.02em;")
        self.lbl_badge.setText(badge_text)
        self.lbl_badge.setStyleSheet(
            f"color: {badge_color}; background: {badge_bg}; "
            f"border: 1px solid {badge_border}; "
            f"border-radius: 5px; padding: 0 10px; font-weight: bold; letter-spacing: 0.1em;"
        )
        self.progress.setStyleSheet(f"""
            QProgressBar {{
                background: rgba(20, 22, 36, 0.9);
                border-radius: 5px;
                border: none;
            }}
            QProgressBar::chunk {{
                border-radius: 5px;
                background: {bar_color};
            }}
        """)

        # Infos secondaires
        info_color = "#5a5c80"
        seuil_txt = f" (seuil : {threshold}%)" if threshold else ""
        self.lbl_nodes.setText(f"{n_nodes} nœud{'s' if n_nodes != 1 else ''} MRD{seuil_txt}")
        self.lbl_nodes.setStyleSheet(
            f"color: {info_color}; background: transparent; font-size: 10px; font-weight: 600;"
        )
        self.lbl_cells.setText(f"{n_cells:,} cellule{'s' if n_cells != 1 else ''}")
        self.lbl_cells.setStyleSheet(
            f"color: {info_color}; background: transparent; font-size: 10px; font-weight: 600;"
        )

    def _set_empty(self) -> None:
        """Affichage vide avant le premier résultat."""
        self.lbl_pct.setText("—")
        self.lbl_badge.setText("EN ATTENTE")
        self.lbl_badge.setStyleSheet(
            f"color: {_SUBTEXT}; background: {_SURFACE1}; "
            f"border-radius: 6px; padding: 0 8px;"
        )
        self.progress.setValue(0)
        self.lbl_nodes.setText("— nœuds")
        self.lbl_cells.setText("— cellules")
