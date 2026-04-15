# -*- coding: utf-8 -*-
"""
log_console.py — Terminal de logs avec coloration syntaxique.

LogConsole est un QPlainTextEdit enrichi avec :
  - Police monospace (Cascadia Code / Fira Code / Consolas)
  - Fond #0a0a14 (terminal sombre)
  - Coloration syntaxique : [INFO] vert, [WARNING] orange, [ERROR] rouge,
    [SUCCESS] cyan, timestamps en gris
  - Méthode append_log(msg) qui parse et colore automatiquement

Design System "Deep Medical Clarity" :
  - [INFO]    → #2ECC71 (Vert Santé)
  - [WARNING] → #F39C12 (Orange Alerte)
  - [ERROR]   → #E74C3C (Rouge Alerte)
  - [SUCCESS] → #00A3FF (Bleu Technologie)
  - Timestamp → #45475a (Gris surface)
"""

from __future__ import annotations

import re
from typing import Optional

from PyQt5.QtWidgets import QPlainTextEdit, QWidget
from PyQt5.QtCore import Qt
from PyQt5.QtGui import (
    QFont,
    QTextCharFormat,
    QColor,
    QTextCursor,
    QBrush,
)

# ── Patterns de coloration ────────────────────────────────────────────

_LOG_RULES: list[tuple[re.Pattern, str]] = [
    # Timestamps ISO : 2024-01-15 12:34:56
    (re.compile(r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}"), "#45475a"),
    # Niveaux de log
    (re.compile(r"\[INFO\]|\[info\]"),       "#2ECC71"),   # Vert santé
    (re.compile(r"\[SUCCESS\]|\[success\]"), "#00A3FF"),   # Bleu technologie
    (re.compile(r"\[WARNING\]|\[warning\]|\[WARN\]|\[warn\]"), "#F39C12"),  # Orange
    (re.compile(r"\[ERROR\]|\[error\]|\[ERR\]|\[err\]"),       "#E74C3C"),  # Rouge
    (re.compile(r"\[DEBUG\]|\[debug\]"),     "#585b70"),   # Gris discret
    # Valeurs numériques importantes (pourcentages, MRD)
    (re.compile(r"\b\d+\.?\d*\s*%"),         "#cba6f7"),   # Mauve — pourcentages
    # Chemins de fichiers
    (re.compile(r"[A-Za-z]:\\[^\s]+|/[^\s]+\.[a-z]+"), "#f9e2af"),  # Jaune doux
    # Étapes pipelines (ex: "Étape 1/5" ou "Step 1 of 5")
    (re.compile(r"[ÉEé]tape\s+\d+\s*/\s*\d+|Step\s+\d+\s+of\s+\d+"), "#89b4fa"),
]


class LogConsole(QPlainTextEdit):
    """
    Terminal de logs avec coloration syntaxique intégrée.

    Usage :
        console = LogConsole()
        console.append_log("[INFO] Pipeline démarré.")
        console.append_log("[ERROR] Fichier introuvable : data.fcs")
    """

    # Nombre max de blocs (lignes) pour éviter une croissance infinie
    MAX_BLOCKS = 5000

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)

        self.setObjectName("logConsole")
        self.setReadOnly(True)
        self.setMaximumBlockCount(self.MAX_BLOCKS)
        self.setPlaceholderText("Les logs du pipeline apparaîtront ici…")

        # Police monospace
        font = QFont()
        font.setFamilies(["Cascadia Code", "Fira Code", "Consolas", "Courier New"])
        font.setPointSize(9)
        font.setStyleHint(QFont.Monospace)
        self.setFont(font)

        # Style de base (complété par QSS global)
        self.setStyleSheet("""
            QPlainTextEdit#logConsole {
                background: #0a0a14;
                border: 1px solid rgba(137, 180, 250, 0.10);
                border-radius: 10px;
                color: #7ec87a;
                padding: 12px;
                selection-background-color: rgba(137, 180, 250, 0.28);
                line-height: 1.5;
            }
            QPlainTextEdit#logConsole QScrollBar:vertical {
                background: #0d0d18;
                width: 8px;
                border: none;
            }
            QPlainTextEdit#logConsole QScrollBar::handle:vertical {
                background: rgba(137, 180, 250, 0.2);
                border-radius: 4px;
                min-height: 20px;
            }
            QPlainTextEdit#logConsole QScrollBar::handle:vertical:hover {
                background: rgba(137, 180, 250, 0.4);
            }
        """)

        # Format de texte par défaut (texte brut)
        self._default_format = QTextCharFormat()
        self._default_format.setForeground(QBrush(QColor("#7ec87a")))

    # ── API publique ──────────────────────────────────────────────────

    def append_log(self, message: str) -> None:
        """
        Ajoute une ligne de log avec coloration syntaxique automatique.
        Scinde les spans par regex et insère chaque segment coloré.
        """
        cursor = self.textCursor()
        cursor.movePosition(QTextCursor.End)

        # Saut de ligne si le document n'est pas vide
        if not self.document().isEmpty():
            cursor.insertText("\n", self._default_format)

        self._insert_colored(cursor, message)

        # Auto-scroll vers le bas
        self.setTextCursor(cursor)
        self.ensureCursorVisible()

    def append(self, text: str) -> None:
        """
        Surcharge de QPlainTextEdit.appendPlainText pour compatibilité avec
        l'API QTextEdit.append() utilisée dans main_window.py.
        Redirige vers append_log() pour appliquer la coloration syntaxique.
        """
        # Retirer les balises HTML simples si le texte vient d'un insertHtml
        clean = re.sub(r"<[^>]+>", "", text)
        self.append_log(clean)

    def clear_logs(self) -> None:
        """Efface tous les logs."""
        self.clear()

    # ── Coloration syntaxique ─────────────────────────────────────────

    def _insert_colored(self, cursor: QTextCursor, line: str) -> None:
        """
        Insère une ligne de log en colorant les tokens correspondant aux règles.
        Algorithme : on calcule tous les spans qui matchent, on les trie,
        on insère segment par segment.
        """
        if not line:
            cursor.insertText("", self._default_format)
            return

        # Collecter tous les spans (start, end, color)
        spans: list[tuple[int, int, str]] = []
        for pattern, color in _LOG_RULES:
            for m in pattern.finditer(line):
                spans.append((m.start(), m.end(), color))

        if not spans:
            # Ligne sans token spécial : couleur de base selon niveau
            fmt = self._format_for_line(line)
            cursor.insertText(line, fmt)
            return

        # Trier et dédupliquer (priorité au premier match en cas de chevauchement)
        spans.sort(key=lambda s: s[0])
        merged: list[tuple[int, int, str]] = []
        last_end = 0
        for start, end, color in spans:
            if start >= last_end:
                merged.append((start, end, color))
                last_end = end

        # Couleur de fond de ligne (selon niveau log)
        base_fmt = self._format_for_line(line)

        pos = 0
        for start, end, color in merged:
            # Segment avant le match
            if pos < start:
                cursor.insertText(line[pos:start], base_fmt)
            # Segment coloré
            fmt = QTextCharFormat()
            fmt.setForeground(QBrush(QColor(color)))
            fmt.setFont(self.font())
            cursor.insertText(line[start:end], fmt)
            pos = end

        # Reste après le dernier match
        if pos < len(line):
            cursor.insertText(line[pos:], base_fmt)

    def _format_for_line(self, line: str) -> QTextCharFormat:
        """Détermine la couleur de base d'une ligne selon son niveau de log."""
        fmt = QTextCharFormat()
        fmt.setFont(self.font())
        line_up = line.upper()
        if "[ERROR]" in line_up or "[ERR]" in line_up:
            fmt.setForeground(QBrush(QColor("#E74C3C")))
        elif "[WARNING]" in line_up or "[WARN]" in line_up:
            fmt.setForeground(QBrush(QColor("#F39C12")))
        elif "[SUCCESS]" in line_up:
            fmt.setForeground(QBrush(QColor("#00A3FF")))
        elif "[DEBUG]" in line_up:
            fmt.setForeground(QBrush(QColor("#585b70")))
        else:
            fmt.setForeground(QBrush(QColor("#7ec87a")))
        return fmt
