# -*- coding: utf-8 -*-
"""
mrd_node_table.py — Tableau haute-performance des nœuds SOM MRD positifs.

Architecture : QTableView + MRDNodeTableModel (QAbstractTableModel)
  - Empreinte mémoire quasi-nulle : le modèle ne copie pas les données.
  - Tri natif via QSortFilterProxyModel.
  - Filtre par méthode via le même proxy.
  - Aucun QTableWidgetItem alloué : les cellules sont rendues à la volée
    par data() / headerData().

Colonnes :
    Nœud SOM | Méthodes | % Sain (nœud) | % Patho (nœud) | Cellules patho | Total nœud
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from PyQt5.QtCore import (
    Qt,
    QAbstractTableModel,
    QModelIndex,
    QSortFilterProxyModel,
    QVariant,
)
from PyQt5.QtGui import QColor, QFont, QBrush
from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QTableView,
    QHeaderView,
    QComboBox,
    QAbstractItemView,
    QSizePolicy,
)


# ── Palette locale (identique à styles.COLORS) ──────────────────────────────
_C = {
    "surface0": "#313244",
    "surface1": "#45475a",
    "base":     "#1e1e2e",
    "crust":    "#11111b",
    "text":     "#cdd6f4",
    "subtext":  "#a6adc8",
    "red":      "#f38ba8",
    "green":    "#a6e3a1",
    "yellow":   "#f9e2af",
    "blue":     "#89b4fa",
    "mauve":    "#cba6f7",
    "overlay0": "#6c7086",
}

# Colonnes
_COL_NODE    = 0
_COL_METHODS = 1
_COL_PCT_SAIN  = 2
_COL_PCT_PATHO = 3
_COL_N_PATHO   = 4
_COL_TOTAL     = 5

_HEADERS = [
    "Nœud SOM",
    "Méthodes",
    "% Sain",
    "% Patho",
    "Cellules patho",
    "Total nœud",
]

# Méthodes → clé booléenne dans le dict node
_METHOD_FLAG: Dict[str, str] = {
    "JF":  "is_mrd_jf",
    "Flo": "is_mrd_flo",
    "ELN": "is_mrd_eln",
}


# ════════════════════════════════════════════════════════════════════
# Modèle de données
# ════════════════════════════════════════════════════════════════════


class MRDNodeTableModel(QAbstractTableModel):
    """
    Modèle léger wrappant une liste de dicts de nœuds MRD.

    Aucune copie n'est faite : data() lit directement dans self._nodes.
    Le tri et le filtrage sont délégués à QSortFilterProxyModel.
    """

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._nodes: List[Dict[str, Any]] = []
        self._font_bold = QFont("Segoe UI", 9, QFont.Bold)
        self._font_mono = QFont("Cascadia Code", 9)

    # ── Interface publique ───────────────────────────────────────────

    def load(self, nodes: List[Dict[str, Any]]) -> None:
        """Remplace les données et notifie la vue."""
        self.beginResetModel()
        self._nodes = nodes
        self.endResetModel()

    def node_at(self, row: int) -> Optional[Dict[str, Any]]:
        if 0 <= row < len(self._nodes):
            return self._nodes[row]
        return None

    # ── QAbstractTableModel interface ───────────────────────────────

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:
        if parent.isValid():
            return 0
        return len(self._nodes)

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return len(_HEADERS)

    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.DisplayRole,
    ) -> Any:
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return _HEADERS[section]
        if orientation == Qt.Horizontal and role == Qt.ForegroundRole:
            return QBrush(QColor(_C["overlay0"]))
        return QVariant()

    def data(self, index: QModelIndex, role: int = Qt.DisplayRole) -> Any:
        if not index.isValid():
            return QVariant()

        row, col = index.row(), index.column()
        if row >= len(self._nodes):
            return QVariant()

        node = self._nodes[row]

        # ── DisplayRole ─────────────────────────────────────────────
        if role == Qt.DisplayRole:
            if col == _COL_NODE:
                return str(node.get("node_id", "?"))
            if col == _COL_METHODS:
                parts = [m for m, key in _METHOD_FLAG.items() if node.get(key)]
                return " · ".join(parts) if parts else "—"
            if col == _COL_PCT_SAIN:
                return f"{node.get('pct_sain', 0.0):.1f} %"
            if col == _COL_PCT_PATHO:
                return f"{node.get('pct_patho', 0.0):.1f} %"
            if col == _COL_N_PATHO:
                n = node.get("n_patho", 0)
                return f"{n:,}"
            if col == _COL_TOTAL:
                n = node.get("n_cells", 0)
                return f"{n:,}"

        # ── TextAlignmentRole ────────────────────────────────────────
        if role == Qt.TextAlignmentRole:
            if col == _COL_METHODS:
                return Qt.AlignCenter | Qt.AlignVCenter
            return Qt.AlignCenter | Qt.AlignVCenter

        # ── ForegroundRole ───────────────────────────────────────────
        if role == Qt.ForegroundRole:
            if col == _COL_METHODS:
                return QBrush(QColor(_C["red"]))
            if col == _COL_PCT_SAIN:
                return QBrush(QColor(_C["green"]))
            if col == _COL_PCT_PATHO:
                return QBrush(QColor(_C["red"]))
            if col in (_COL_N_PATHO, _COL_TOTAL):
                return QBrush(QColor(_C["subtext"]))
            return QBrush(QColor(_C["text"]))

        # ── BackgroundRole (ligne alternée) ──────────────────────────
        if role == Qt.BackgroundRole:
            if row % 2 == 0:
                return QBrush(QColor(14, 15, 24, 220))
            return QBrush(QColor(24, 26, 42, 180))

        # ── FontRole ─────────────────────────────────────────────────
        if role == Qt.FontRole:
            if col == _COL_NODE:
                return self._font_bold
            if col in (_COL_N_PATHO, _COL_TOTAL):
                return self._font_mono
            return QVariant()

        # ── SortRole (tri numérique) ─────────────────────────────────
        if role == Qt.UserRole:
            if col == _COL_NODE:
                return int(node.get("node_id", 0))
            if col == _COL_PCT_SAIN:
                return float(node.get("pct_sain", 0.0))
            if col == _COL_PCT_PATHO:
                return float(node.get("pct_patho", 0.0))
            if col == _COL_N_PATHO:
                return int(node.get("n_patho", 0))
            if col == _COL_TOTAL:
                return int(node.get("n_cells", 0))
            return QVariant()

        return QVariant()

    def sort(self, column: int, order: Qt.SortOrder = Qt.AscendingOrder) -> None:
        """Tri in-place pour support de QTableView.sortByColumn()."""
        reverse = order == Qt.DescendingOrder

        def key_fn(node: Dict[str, Any]) -> Any:
            if column == _COL_NODE:
                return int(node.get("node_id", 0))
            if column == _COL_METHODS:
                parts = [m for m, k in _METHOD_FLAG.items() if node.get(k)]
                return " ".join(parts)
            if column == _COL_PCT_SAIN:
                return float(node.get("pct_sain", 0.0))
            if column == _COL_PCT_PATHO:
                return float(node.get("pct_patho", 0.0))
            if column == _COL_N_PATHO:
                return int(node.get("n_patho", 0))
            if column == _COL_TOTAL:
                return int(node.get("n_cells", 0))
            return 0

        self.beginResetModel()
        self._nodes.sort(key=key_fn, reverse=reverse)
        self.endResetModel()


# ════════════════════════════════════════════════════════════════════
# Proxy de filtre par méthode
# ════════════════════════════════════════════════════════════════════


class MRDMethodFilterProxy(QSortFilterProxyModel):
    """
    Filtre les lignes selon la méthode sélectionnée.
    method="" → toutes les lignes.
    """

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._method: str = ""

    def set_method(self, method: str) -> None:
        self._method = method
        self.invalidateFilter()

    def filterAcceptsRow(
        self, source_row: int, source_parent: QModelIndex
    ) -> bool:
        if not self._method:
            return True
        model = self.sourceModel()
        if not isinstance(model, MRDNodeTableModel):
            return True
        node = model.node_at(source_row)
        if node is None:
            return False
        flag_key = _METHOD_FLAG.get(self._method, "")
        if not flag_key:
            return True
        return bool(node.get(flag_key, False))


# ════════════════════════════════════════════════════════════════════
# Widget composite
# ════════════════════════════════════════════════════════════════════


class MRDNodeTable(QWidget):
    """
    Widget composite : en-tête + filtre méthode + QTableView haute-performance.

    Interface publique :
        load_nodes(nodes, available_methods)  → charge les données
        clear()                               → vide tout
    """

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._build_ui()

    # ── Construction ────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(8)

        # ── En-tête ──────────────────────────────────────────────────
        header_widget = QWidget()
        header_widget.setStyleSheet("""
            QWidget {
                background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                    stop:0 rgba(30, 32, 52, 0.8), stop:1 rgba(22, 24, 40, 0.8));
                border-radius: 8px 8px 0px 0px;
                border: 1px solid rgba(137, 180, 250, 0.12);
                border-bottom: none;
            }
        """)
        header = QHBoxLayout(header_widget)
        header.setContentsMargins(14, 10, 14, 10)
        header.setSpacing(10)

        lbl = QLabel("NŒUDS SOM MRD POSITIFS")
        lbl.setObjectName("sectionLabel")
        lbl.setStyleSheet(
            "color: #5a5c80; font-size: 9pt; font-weight: 700; "
            "letter-spacing: 0.1em; background: transparent;"
        )
        header.addWidget(lbl)
        header.addStretch()

        lbl_filter = QLabel("Filtre :")
        lbl_filter.setStyleSheet(
            f"color: {_C['overlay0']}; font-size: 9pt; background: transparent;"
        )
        header.addWidget(lbl_filter)

        self.combo_filter = QComboBox()
        self.combo_filter.setMinimumWidth(140)
        self.combo_filter.setMaximumWidth(180)
        self.combo_filter.currentTextChanged.connect(self._on_filter_changed)
        header.addWidget(self.combo_filter)

        root.addWidget(header_widget)

        # ── Modèle + Proxy ───────────────────────────────────────────
        self._model = MRDNodeTableModel(self)
        self._proxy = MRDMethodFilterProxy(self)
        self._proxy.setSourceModel(self._model)
        self._proxy.setSortRole(Qt.UserRole)

        # ── Vue ──────────────────────────────────────────────────────
        self._view = QTableView()
        self._view.setModel(self._proxy)
        self._view.setSortingEnabled(True)
        self._view.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self._view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._view.setSelectionMode(QAbstractItemView.SingleSelection)
        self._view.setAlternatingRowColors(False)   # géré par BackgroundRole
        self._view.verticalHeader().setVisible(False)
        self._view.setShowGrid(False)
        self._view.setWordWrap(False)

        # Hauteur des rangées
        self._view.verticalHeader().setDefaultSectionSize(38)

        # Colonnes
        hh = self._view.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.Interactive)
        hh.setSectionResizeMode(_COL_NODE,      QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(_COL_METHODS,   QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(_COL_PCT_SAIN,  QHeaderView.Stretch)
        hh.setSectionResizeMode(_COL_PCT_PATHO, QHeaderView.Stretch)
        hh.setSectionResizeMode(_COL_N_PATHO,   QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(_COL_TOTAL,     QHeaderView.ResizeToContents)
        hh.setHighlightSections(False)
        hh.setStretchLastSection(False)

        # Tri par défaut : % Patho décroissant
        self._view.sortByColumn(_COL_PCT_PATHO, Qt.DescendingOrder)

        root.addWidget(self._view, 1)

        # ── Label état vide ──────────────────────────────────────────
        self._lbl_empty = QLabel("Aucun nœud MRD positif détecté pour cette sélection")
        self._lbl_empty.setAlignment(Qt.AlignCenter)
        self._lbl_empty.setStyleSheet(
            f"color: {_C['overlay0']}; font-style: italic; padding: 24px; font-size: 10pt;"
        )
        self._lbl_empty.hide()
        root.addWidget(self._lbl_empty)

        # Compteur de lignes
        self._lbl_count = QLabel("")
        self._lbl_count.setAlignment(Qt.AlignRight)
        self._lbl_count.setStyleSheet(
            f"color: {_C['overlay0']}; font-size: 8pt; padding: 2px 4px;"
        )
        root.addWidget(self._lbl_count)

        # Mise à jour du compteur à chaque changement du proxy
        self._proxy.rowsInserted.connect(self._refresh_empty_state)
        self._proxy.rowsRemoved.connect(self._refresh_empty_state)
        self._proxy.modelReset.connect(self._refresh_empty_state)
        self._model.modelReset.connect(self._refresh_empty_state)

    # ── Interface publique ───────────────────────────────────────────

    def load_nodes(
        self,
        nodes: List[Dict[str, Any]],
        available_methods: List[str],
    ) -> None:
        """Charge les nœuds MRD, peuple le filtre, déclenche l'affichage."""
        # Filtre
        self.combo_filter.blockSignals(True)
        self.combo_filter.clear()
        self.combo_filter.addItem("Toutes les méthodes")
        for m in available_methods:
            self.combo_filter.addItem(m)
        self.combo_filter.setCurrentIndex(0)
        self.combo_filter.blockSignals(False)

        # Données
        self._proxy.set_method("")
        self._model.load(nodes)

        # Tri par défaut
        self._view.sortByColumn(_COL_PCT_PATHO, Qt.DescendingOrder)
        self._refresh_empty_state()

    def clear(self) -> None:
        self._model.load([])
        self.combo_filter.blockSignals(True)
        self.combo_filter.clear()
        self.combo_filter.blockSignals(False)
        self._proxy.set_method("")
        self._refresh_empty_state()

    # ── Slots internes ───────────────────────────────────────────────

    def _on_filter_changed(self, text: str) -> None:
        if text in ("", "Toutes les méthodes"):
            self._proxy.set_method("")
        else:
            # Normalise "ELN 2025" → "ELN"
            clean = text.split()[0] if text else ""
            self._proxy.set_method(clean)
        self._refresh_empty_state()

    def _refresh_empty_state(self) -> None:
        visible = self._proxy.rowCount()
        total = self._model.rowCount()

        if visible == 0:
            self._view.hide()
            self._lbl_empty.show()
            self._lbl_count.setText("")
        else:
            self._lbl_empty.hide()
            self._view.show()
            if visible == total:
                self._lbl_count.setText(f"{visible} nœud(s) MRD")
            else:
                self._lbl_count.setText(f"{visible}/{total} nœud(s) affiché(s)")
