# -*- coding: utf-8 -*-
"""
settings_card.py — Carte de paramétrage "Deep Medical Clarity".

SettingsCard est un QFrame avec :
  - Un titre stylisé en haut
  - Un fond légèrement plus clair que le fond principal
  - Des coins arrondis (12px)
  - Une bordure fine semi-transparente

Usage :
    card = SettingsCard("Paramètres FlowSOM", icon_name="fa5s.sliders-h")
    card.body_layout.addWidget(mon_widget)
"""

from __future__ import annotations

from typing import Optional

from PyQt5.QtWidgets import (
    QFrame,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QWidget,
    QSizePolicy,
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont

try:
    import qtawesome as qta
    _QTA = True
except ImportError:
    _QTA = False


class SettingsCard(QFrame):
    """
    Carte de regroupement de paramètres avec titre, icône et corps scrollable.

    Attributes :
        body_layout (QVBoxLayout) : layout dans lequel ajouter les widgets enfants
    """

    def __init__(
        self,
        title: str,
        icon_name: str = "",
        subtitle: str = "",
        accent: str = "#00A3FF",
        parent: Optional[QWidget] = None,
    ) -> None:
        super().__init__(parent)

        self._accent = accent
        self.setObjectName("settingsCard")
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)

        # Style inline pour les variantes d'accent (override QSS)
        self.setStyleSheet(f"""
            QFrame#settingsCard {{
                background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                    stop:0 rgba(35, 37, 55, 0.95),
                    stop:1 rgba(27, 29, 44, 0.90));
                border: 1px solid rgba(137, 180, 250, 0.12);
                border-top: 1px solid rgba({self._hex_to_rgb(accent)}, 0.25);
                border-radius: 12px;
            }}
        """)

        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ── En-tête de la carte ───────────────────────────────────────
        header = QWidget(self)
        header.setObjectName("cardHeader")
        header.setStyleSheet(f"""
            QWidget#cardHeader {{
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0 rgba({self._hex_to_rgb(accent)}, 0.12),
                    stop:1 transparent);
                border-bottom: 1px solid rgba(137, 180, 250, 0.08);
                border-top-left-radius: 12px;
                border-top-right-radius: 12px;
            }}
        """)
        h_layout = QHBoxLayout(header)
        h_layout.setContentsMargins(16, 12, 16, 12)
        h_layout.setSpacing(10)

        # Icône optionnelle
        if icon_name and _QTA:
            try:
                icon_lbl = QLabel(header)
                ico = qta.icon(icon_name, color=accent)
                from PyQt5.QtCore import QSize
                icon_lbl.setPixmap(ico.pixmap(QSize(16, 16)))
                icon_lbl.setFixedSize(20, 20)
                icon_lbl.setStyleSheet("background: transparent;")
                h_layout.addWidget(icon_lbl)
            except Exception:
                pass

        # Titre
        title_lbl = QLabel(title, header)
        title_lbl.setFont(QFont("Segoe UI", 10, QFont.Bold))
        title_lbl.setStyleSheet(
            f"color: {accent}; background: transparent; letter-spacing: 0.02em;"
        )
        h_layout.addWidget(title_lbl)

        # Sous-titre optionnel
        if subtitle:
            sub_lbl = QLabel(f" — {subtitle}", header)
            sub_lbl.setStyleSheet(
                "color: #585b70; font-size: 8.5pt; background: transparent;"
            )
            h_layout.addWidget(sub_lbl)

        h_layout.addStretch()
        root.addWidget(header)

        # ── Corps de la carte ─────────────────────────────────────────
        body = QWidget(self)
        body.setObjectName("cardBody")
        body.setStyleSheet("QWidget#cardBody { background: transparent; }")

        self.body_layout = QVBoxLayout(body)
        self.body_layout.setContentsMargins(16, 14, 16, 16)
        self.body_layout.setSpacing(10)
        root.addWidget(body)

    # ── Helpers ───────────────────────────────────────────────────────

    @staticmethod
    def _hex_to_rgb(hex_color: str) -> str:
        """Convertit #RRGGBB → 'R, G, B' pour usage dans rgba()."""
        h = hex_color.lstrip("#")
        if len(h) == 6:
            r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
            return f"{r}, {g}, {b}"
        return "137, 180, 250"  # fallback bleu Catppuccin
