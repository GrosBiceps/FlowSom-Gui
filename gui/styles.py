# -*- coding: utf-8 -*-
"""
styles.py — Thème Catppuccin Mocha / Flat-Design Enterprise pour FlowSomAnalyzerPro.

Palette Catppuccin Mocha :
    Base       #1e1e2e    Mantle   #181825    Crust    #11111b
    Surface0   #313244    Surface1 #45475a    Surface2 #585b70
    Overlay0   #6c7086    Text     #cdd6f4    Subtext  #a6adc8
    Blue       #89b4fa    Lavender #b4befe    Mauve    #cba6f7
    Pink       #f5c2e7    Red      #f38ba8    Peach    #fab387
    Green      #a6e3a1    Teal     #94e2d5    Sky      #89dceb
    Sapphire   #74c7ec

Design tokens :
  - Font principale  : Segoe UI → Inter → Roboto (ordre de préférence Windows)
  - Font monospace   : Cascadia Code → Consolas → Fira Code
  - Radius standard  : 8px (inputs), 10px (cards/groupbox), 6px (badges)
  - Transition        : simulé via hover contrast (QSS ne supporte pas CSS transitions)
"""

# ---------------------------------------------------------------------------
# Couleurs réutilisables
# ---------------------------------------------------------------------------
COLORS = {
    "base": "#1e1e2e",
    "mantle": "#181825",
    "crust": "#11111b",
    "surface0": "#313244",
    "surface1": "#45475a",
    "surface2": "#585b70",
    "overlay0": "#6c7086",
    "text": "#cdd6f4",
    "subtext": "#a6adc8",
    "blue": "#89b4fa",
    "lavender": "#b4befe",
    "mauve": "#cba6f7",
    "pink": "#f5c2e7",
    "red": "#f38ba8",
    "peach": "#fab387",
    "green": "#a6e3a1",
    "teal": "#94e2d5",
    "sky": "#89dceb",
    "sapphire": "#74c7ec",
    # Accents additionnels
    "blue_bright": "#a8c8ff",
    "green_bright": "#b8f0b4",
    "red_bright": "#f7a8bf",
}

# ---------------------------------------------------------------------------
# Indicateurs de progression du Stepper (couleurs sémantiques)
# ---------------------------------------------------------------------------
STEP_COLORS = {
    "done": "#a6e3a1",  # vert — étape terminée
    "active": "#89b4fa",  # bleu — étape active
    "pending": "#45475a",  # gris — étape future
    "error": "#f38ba8",  # rouge — erreur
}

# ---------------------------------------------------------------------------
# QSS global — Flat Design Enterprise / Catppuccin Mocha
# ---------------------------------------------------------------------------
STYLESHEET = """
/* ════════════════════════════════════════════════════════════════════
   BASE
   ════════════════════════════════════════════════════════════════════ */

QMainWindow {
    background: #0d0d17;
}

QWidget {
    background-color: transparent;
    color: #cdd6f4;
    font-family: 'Segoe UI', 'Inter', 'Roboto', Arial, sans-serif;
    font-size: 10pt;
}

QWidget#centralWidget {
    background: #0d0d17;
}

/* ── Wizard / Stepper Shell ──────────────────────────────────────── */

QWidget#wizardShell {
    background: #0d0d17;
}

/* Barre latérale des étapes (sidebar) */
QWidget#stepSidebar {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #181825, stop:1 #11111b);
    border-right: 1px solid rgba(137, 180, 250, 0.12);
}

/* Contenu principal (droite du wizard) */
QWidget#stepContent {
    background: #13131f;
    border-radius: 0px;
}

/* Panneau gauche contrôles (mode non-wizard) */
QWidget#leftPanel {
    background: #181825;
    border-right: 1px solid rgba(137, 180, 250, 0.12);
}

QWidget#rightPanel {
    background: #13131f;
}

/* ════════════════════════════════════════════════════════════════════
   STEPPER — étapes de navigation
   ════════════════════════════════════════════════════════════════════ */

QLabel#stepLabel {
    font-size: 9pt;
    font-weight: 600;
    letter-spacing: 0.05em;
    text-transform: uppercase;
    padding: 3px 0px;
}

QLabel#stepTitle {
    font-size: 11pt;
    font-weight: 700;
    color: #cdd6f4;
}

QLabel#stepDone   { color: #a6e3a1; }
QLabel#stepActive { color: #89b4fa; }
QLabel#stepPending{ color: #3a3c52; }

/* Bouton d'étape dans la sidebar */
QPushButton#stepBtn {
    background: transparent;
    border: none;
    border-left: 3px solid transparent;
    border-radius: 0px;
    padding: 11px 16px;
    text-align: left;
    font-size: 10pt;
    font-weight: 500;
    color: #4a4c65;
}

QPushButton#stepBtn:hover {
    background: rgba(137, 180, 250, 0.07);
    color: #8a8cb0;
    border-left-color: rgba(137, 180, 250, 0.25);
}

QPushButton#stepBtnActive {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 rgba(137, 180, 250, 0.18), stop:1 transparent);
    border: none;
    border-left: 3px solid #89b4fa;
    border-radius: 0px;
    padding: 11px 16px;
    text-align: left;
    font-size: 10pt;
    font-weight: 700;
    color: #c8d8fd;
}

QPushButton#stepBtnDone {
    background: transparent;
    border: none;
    border-left: 3px solid rgba(166, 227, 161, 0.6);
    border-radius: 0px;
    padding: 11px 16px;
    text-align: left;
    font-size: 10pt;
    font-weight: 500;
    color: #7ec87a;
}

QPushButton#stepBtnDone:hover {
    background: rgba(166, 227, 161, 0.06);
    color: #a6e3a1;
}

/* ════════════════════════════════════════════════════════════════════
   GROUPBOX — cartes glassmorphism
   ════════════════════════════════════════════════════════════════════ */

QGroupBox {
    font-weight: 700;
    font-size: 10pt;
    letter-spacing: 0.03em;
    border: 1px solid rgba(137, 180, 250, 0.14);
    border-top: 1px solid rgba(137, 180, 250, 0.22);
    border-radius: 12px;
    margin-top: 20px;
    padding: 20px 16px 16px 16px;
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 rgba(55, 56, 76, 0.55), stop:1 rgba(40, 41, 58, 0.45));
}

QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 16px;
    top: 5px;
    padding: 4px 14px;
    color: #a8c8ff;
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 rgba(137, 180, 250, 0.22), stop:1 rgba(137, 180, 250, 0.08));
    border-radius: 6px;
    border: 1px solid rgba(137, 180, 250, 0.2);
    font-size: 10pt;
    font-weight: 700;
}

/* ════════════════════════════════════════════════════════════════════
   DRAG & DROP ZONE (Step 1 — Import)
   ════════════════════════════════════════════════════════════════════ */

QLabel#dropZone {
    background: rgba(137, 180, 250, 0.04);
    border: 2px dashed rgba(137, 180, 250, 0.28);
    border-radius: 14px;
    color: #4a4c65;
    font-size: 10pt;
    padding: 36px 20px;
    qproperty-alignment: AlignCenter;
}

QLabel#dropZone[dragOver="true"] {
    background: rgba(137, 180, 250, 0.11);
    border-color: #89b4fa;
    color: #cdd6f4;
}

QLabel#dropZoneOk {
    background: rgba(166, 227, 161, 0.07);
    border: 2px solid rgba(166, 227, 161, 0.38);
    border-radius: 14px;
    color: #a6e3a1;
    font-size: 10pt;
    padding: 36px 20px;
    qproperty-alignment: AlignCenter;
}

/* ════════════════════════════════════════════════════════════════════
   BOUTONS — hiérarchie claire avec depth
   ════════════════════════════════════════════════════════════════════ */

/* Base (secondaire) */
QPushButton {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 rgba(78, 80, 102, 0.75), stop:1 rgba(58, 60, 80, 0.75));
    border: 1px solid rgba(137, 180, 250, 0.2);
    border-bottom: 1px solid rgba(0, 0, 0, 0.3);
    border-radius: 8px;
    padding: 8px 18px;
    color: #c5cef0;
    font-weight: 600;
    font-size: 10pt;
    min-height: 16px;
}

QPushButton:hover {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 rgba(95, 97, 122, 0.9), stop:1 rgba(72, 74, 96, 0.9));
    border-color: rgba(137, 180, 250, 0.45);
    color: #dde6ff;
}

QPushButton:pressed {
    background: rgba(44, 46, 62, 0.95);
    border-color: rgba(137, 180, 250, 0.5);
    padding-top: 9px;
    padding-bottom: 7px;
}

QPushButton:disabled {
    background: rgba(35, 37, 52, 0.5);
    color: #3a3c52;
    border-color: rgba(69, 71, 90, 0.2);
}

/* Primaire (CTA) — bleu lumineux */
QPushButton#primaryBtn {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #9fc4ff, stop:1 #7aacf7);
    border: 1px solid rgba(137, 180, 250, 0.5);
    border-bottom: 1px solid rgba(80, 120, 210, 0.6);
    border-radius: 8px;
    color: #eef4ff;
    font-weight: 700;
    font-size: 10pt;
    padding: 10px 24px;
}

QPushButton#primaryBtn:hover {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #b4d0ff, stop:1 #8fbbfe);
    border-color: rgba(180, 210, 255, 0.7);
    color: #f4f8ff;
}

QPushButton#primaryBtn:pressed {
    background: #6ea0ec;
    padding-top: 11px;
    padding-bottom: 9px;
    color: #eaf2ff;
}

QPushButton#primaryBtn:disabled {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 rgba(128, 160, 214, 0.5), stop:1 rgba(96, 126, 178, 0.5));
    color: #dbe6ff;
    border-color: rgba(137, 180, 250, 0.52);
}

/* Succès — vert */
QPushButton#successBtn {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #b8f0b4, stop:1 #90d48c);
    border: 1px solid rgba(166, 227, 161, 0.5);
    border-bottom: 1px solid rgba(100, 170, 96, 0.5);
    color: #0d1f0d;
    font-weight: 700;
    padding: 8px 18px;
    border-radius: 8px;
}

QPushButton#successBtn:hover {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #c8f8c4, stop:1 #a0e49c);
}

QPushButton#successBtn:pressed {
    background: #7ec47a;
    padding-top: 9px;
    padding-bottom: 7px;
}

QPushButton#successBtn:disabled {
    background: rgba(166, 227, 161, 0.18);
    color: rgba(13, 31, 13, 0.35);
    border-color: rgba(166, 227, 161, 0.12);
}

/* Danger — rouge */
QPushButton#dangerBtn {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #f7a8bf, stop:1 #e87090);
    border: 1px solid rgba(243, 139, 168, 0.5);
    border-bottom: 1px solid rgba(180, 60, 90, 0.5);
    color: #1f0d10;
    font-weight: 700;
    padding: 8px 18px;
    border-radius: 8px;
}

QPushButton#dangerBtn:hover {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #fbbace, stop:1 #f080a4);
}

QPushButton#dangerBtn:pressed {
    background: #d65e80;
    padding-top: 9px;
    padding-bottom: 7px;
}

QPushButton#dangerBtn:disabled {
    background: rgba(243, 139, 168, 0.18);
    color: rgba(31, 13, 16, 0.35);
    border-color: rgba(243, 139, 168, 0.12);
}

/* Export — mauve */
QPushButton#exportBtn {
    background: rgba(203, 166, 247, 0.13);
    border: 1px solid rgba(203, 166, 247, 0.35);
    color: #cba6f7;
    font-weight: 600;
    padding: 8px 16px;
    border-radius: 8px;
}

QPushButton#exportBtn:hover {
    background: rgba(203, 166, 247, 0.22);
    border-color: rgba(203, 166, 247, 0.6);
    color: #dbbeff;
}

QPushButton#exportBtn:pressed {
    background: rgba(203, 166, 247, 0.3);
}

QPushButton#exportBtn:disabled {
    background: rgba(203, 166, 247, 0.05);
    color: rgba(203, 166, 247, 0.25);
    border-color: rgba(203, 166, 247, 0.1);
}

/* Ghost / plat — outline uniquement */
QPushButton#ghostBtn {
    background: transparent;
    border: 1px solid rgba(137, 180, 250, 0.28);
    color: #89b4fa;
    font-weight: 600;
    padding: 7px 16px;
    border-radius: 8px;
}

QPushButton#ghostBtn:hover {
    background: rgba(137, 180, 250, 0.1);
    border-color: rgba(137, 180, 250, 0.55);
    color: #b4d0ff;
}

QPushButton#ghostBtn:pressed {
    background: rgba(137, 180, 250, 0.15);
    border-color: rgba(137, 180, 250, 0.7);
}

/* Icône seule — carré compact */
QPushButton#iconBtn {
    background: rgba(60, 62, 82, 0.55);
    border: 1px solid rgba(137, 180, 250, 0.15);
    border-radius: 8px;
    padding: 6px;
    min-width: 32px;
    max-width: 36px;
    min-height: 32px;
    max-height: 36px;
}

QPushButton#iconBtn:hover {
    background: rgba(137, 180, 250, 0.15);
    border-color: rgba(137, 180, 250, 0.4);
}

/* ════════════════════════════════════════════════════════════════════
   LABELS — hiérarchie typographique
   ════════════════════════════════════════════════════════════════════ */

QLabel {
    color: #cdd6f4;
    background: transparent;
}

QLabel#titleLabel {
    font-size: 20pt;
    font-weight: 700;
    color: #c8d8fd;
    letter-spacing: -0.03em;
}

QLabel#subtitleLabel {
    font-size: 10pt;
    color: #4a4c65;
    font-weight: 400;
}

QLabel#sectionLabel {
    font-size: 9pt;
    font-weight: 700;
    color: #7a7c98;
    letter-spacing: 0.08em;
    text-transform: uppercase;
}

QLabel#cardLabel {
    font-size: 8pt;
    color: #4a4c65;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.1em;
}

QLabel#fileLabel {
    font-size: 9pt;
    color: #a6adc8;
    padding: 6px 12px;
    background: rgba(40, 42, 60, 0.7);
    border: 1px solid rgba(137, 180, 250, 0.15);
    border-radius: 7px;
}

QLabel#statusSuccess {
    color: #a6e3a1;
    font-weight: 600;
}

QLabel#statusError {
    color: #f38ba8;
    font-weight: 600;
}

QLabel#statusWarning {
    color: #f9e2af;
    font-weight: 600;
}

/* Badge pill */
QLabel#badge {
    background: rgba(137, 180, 250, 0.16);
    color: #a8c8ff;
    border-radius: 10px;
    border: 1px solid rgba(137, 180, 250, 0.25);
    padding: 2px 10px;
    font-size: 8pt;
    font-weight: 700;
}

QLabel#badgeGreen {
    background: rgba(166, 227, 161, 0.16);
    color: #a6e3a1;
    border-radius: 10px;
    border: 1px solid rgba(166, 227, 161, 0.25);
    padding: 2px 10px;
    font-size: 8pt;
    font-weight: 700;
}

QLabel#badgeRed {
    background: rgba(243, 139, 168, 0.16);
    color: #f38ba8;
    border-radius: 10px;
    border: 1px solid rgba(243, 139, 168, 0.25);
    padding: 2px 10px;
    font-size: 8pt;
    font-weight: 700;
}

/* ════════════════════════════════════════════════════════════════════
   VALIDATION — bordure rouge si invalide
   ════════════════════════════════════════════════════════════════════ */

QSpinBox[invalid="true"],
QDoubleSpinBox[invalid="true"],
QLineEdit[invalid="true"] {
    border: 1px solid #f38ba8;
    background: rgba(243, 139, 168, 0.07);
}

/* ════════════════════════════════════════════════════════════════════
   INPUTS
   ════════════════════════════════════════════════════════════════════ */

QSpinBox, QDoubleSpinBox, QComboBox, QLineEdit {
    background: rgba(28, 30, 45, 0.85);
    border: 1px solid rgba(137, 180, 250, 0.18);
    border-radius: 8px;
    padding: 7px 10px;
    color: #cdd6f4;
    min-height: 16px;
    font-size: 10pt;
    selection-background-color: rgba(137, 180, 250, 0.35);
    selection-color: #cdd6f4;
}

QSpinBox:hover, QDoubleSpinBox:hover, QComboBox:hover, QLineEdit:hover {
    border-color: rgba(137, 180, 250, 0.35);
    background: rgba(32, 34, 52, 0.95);
}

QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus, QLineEdit:focus {
    border: 1px solid rgba(137, 180, 250, 0.7);
    background: rgba(22, 24, 40, 1.0);
    outline: none;
}

QSpinBox::up-button, QSpinBox::down-button,
QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {
    background: rgba(137, 180, 250, 0.1);
    border: none;
    border-radius: 4px;
    width: 16px;
    margin: 2px;
}

QSpinBox::up-button:hover, QSpinBox::down-button:hover,
QDoubleSpinBox::up-button:hover, QDoubleSpinBox::down-button:hover {
    background: rgba(137, 180, 250, 0.25);
}

QComboBox::drop-down {
    border: none;
    width: 26px;
    background: rgba(137, 180, 250, 0.1);
    border-top-right-radius: 8px;
    border-bottom-right-radius: 8px;
}

/* ── Popup ComboBox (fix fond blanc Windows) ─────────────────────── */

QComboBoxPrivateContainer {
    background: #252637;
    border: 1px solid rgba(137, 180, 250, 0.25);
    margin: 0px;
    padding: 0px;
}

QComboBoxPrivateScroller {
    background: #252637;
    border: none;
    margin: 0px;
    padding: 0px;
    min-height: 0px;
}

QComboBoxPrivateContainer QWidget {
    background: #252637;
    border: none;
}

QComboBox QAbstractItemView {
    background: #252637;
    border: 1px solid #252637;
    selection-background-color: rgba(137, 180, 250, 0.22);
    selection-color: #cdd6f4;
    outline: none;
    padding: 2px;
}

QComboBox QAbstractItemView::item {
    background: #252637;
    color: #cdd6f4;
    padding: 7px 12px;
    border: none;
    min-height: 22px;
    border-radius: 4px;
}

QComboBox QAbstractItemView::item:hover,
QComboBox QAbstractItemView::item:focus {
    background: rgba(137, 180, 250, 0.15);
}

QComboBox QAbstractItemView::item:selected {
    background: rgba(137, 180, 250, 0.28);
    color: #cdd6f4;
}

QComboBox QScrollBar:vertical {
    background: #252637;
    width: 7px;
    border: none;
}

QComboBox QScrollBar::handle:vertical {
    background: rgba(137, 180, 250, 0.3);
    border-radius: 3px;
    min-height: 20px;
}

QComboBox QFrame {
    background: #252637;
    border: none;
}

/* ════════════════════════════════════════════════════════════════════
   CHECKBOX / RADIOBUTTON
   ════════════════════════════════════════════════════════════════════ */

QCheckBox {
    spacing: 10px;
    color: #b8c2e8;
    font-size: 10pt;
}

QCheckBox::indicator {
    width: 18px;
    height: 18px;
    border-radius: 5px;
    border: 1px solid rgba(137, 180, 250, 0.3);
    background: rgba(28, 30, 48, 0.85);
}

QCheckBox::indicator:hover {
    border-color: rgba(137, 180, 250, 0.55);
    background: rgba(35, 38, 60, 0.95);
}

QCheckBox::indicator:checked {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 #9fc4ff, stop:1 #7aacf7);
    border-color: #89b4fa;
    image: none;
}

QCheckBox::indicator:checked:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 #b4d0ff, stop:1 #8fbbfe);
}

QCheckBox:disabled {
    color: #2e3048;
}

QCheckBox::indicator:disabled {
    border-color: rgba(137, 180, 250, 0.08);
    background: rgba(28, 30, 48, 0.35);
}

/* ════════════════════════════════════════════════════════════════════
   PROGRESSBAR — pipeline
   ════════════════════════════════════════════════════════════════════ */

QProgressBar {
    border: none;
    border-radius: 8px;
    background: rgba(28, 30, 48, 0.8);
    text-align: center;
    color: #cdd6f4;
    font-weight: 700;
    font-size: 9pt;
    min-height: 18px;
    max-height: 18px;
}

QProgressBar::chunk {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 #7aacf7, stop:0.45 #a088f0, stop:0.75 #c084c0, stop:1 #90d490);
    border-radius: 8px;
}

/* ════════════════════════════════════════════════════════════════════
   LOG CONSOLE — step Exécution
   ════════════════════════════════════════════════════════════════════ */

QTextEdit#logConsole {
    background: #0a0a14;
    border: 1px solid rgba(137, 180, 250, 0.1);
    border-radius: 10px;
    color: #88e888;
    padding: 12px;
    font-family: 'Cascadia Code', 'Fira Code', 'Consolas', monospace;
    font-size: 9pt;
    selection-background-color: rgba(137, 180, 250, 0.28);
    line-height: 1.5;
}

QTextEdit {
    background: rgba(10, 10, 20, 0.9);
    border: 1px solid rgba(137, 180, 250, 0.12);
    border-radius: 10px;
    color: #88e888;
    padding: 12px;
    font-family: 'Cascadia Code', 'Fira Code', 'Consolas', monospace;
    font-size: 9pt;
    selection-background-color: rgba(137, 180, 250, 0.28);
}

/* ════════════════════════════════════════════════════════════════════
   STATUSBAR
   ════════════════════════════════════════════════════════════════════ */

QStatusBar {
    background: #11111b;
    color: #4a4c65;
    border-top: 1px solid rgba(137, 180, 250, 0.08);
    font-size: 9pt;
    padding: 4px 12px;
}

/* ════════════════════════════════════════════════════════════════════
   TABWIDGET — onglets résultats (panneau droit)
   ════════════════════════════════════════════════════════════════════ */

QTabWidget::pane {
    border: 1px solid rgba(137, 180, 250, 0.14);
    border-top: none;
    border-radius: 0px 0px 12px 12px;
    background: rgba(19, 19, 31, 0.98);
    padding-top: 0px;
}

QTabWidget::tab-bar {
    left: 0px;
}

QTabBar {
    qproperty-drawBase: 0;
    background: transparent;
}

QTabBar::tab {
    background: rgba(30, 32, 48, 0.7);
    border: none;
    border-bottom: 2px solid transparent;
    padding: 10px 18px;
    margin-right: 2px;
    color: #4a4c65;
    font-weight: 500;
    font-size: 9pt;
    min-width: 80px;
}

QTabBar::tab:selected {
    background: rgba(13, 13, 23, 0.95);
    color: #a8c8ff;
    font-weight: 700;
    border-bottom: 2px solid #89b4fa;
}

QTabBar::tab:hover:!selected {
    background: rgba(40, 42, 62, 0.75);
    color: #7a7c98;
    border-bottom: 2px solid rgba(137, 180, 250, 0.2);
}

/* ════════════════════════════════════════════════════════════════════
   SCROLLAREA / SCROLLBAR
   ════════════════════════════════════════════════════════════════════ */

QScrollArea {
    border: none;
    background: transparent;
}

QScrollBar:vertical {
    background: transparent;
    width: 6px;
    border-radius: 3px;
    margin: 4px 2px;
}

QScrollBar::handle:vertical {
    background: rgba(137, 180, 250, 0.18);
    border-radius: 3px;
    min-height: 28px;
}

QScrollBar::handle:vertical:hover {
    background: rgba(137, 180, 250, 0.38);
}

QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0;
}

QScrollBar:horizontal {
    background: transparent;
    height: 6px;
    border-radius: 3px;
    margin: 2px 4px;
}

QScrollBar::handle:horizontal {
    background: rgba(137, 180, 250, 0.18);
    border-radius: 3px;
    min-width: 28px;
}

QScrollBar::handle:horizontal:hover {
    background: rgba(137, 180, 250, 0.38);
}

QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
    width: 0;
}

/* ════════════════════════════════════════════════════════════════════
   TABLE — QTableView + QTableWidget
   ════════════════════════════════════════════════════════════════════ */

QTableView, QTableWidget {
    background: rgba(14, 15, 24, 0.95);
    border: 1px solid rgba(137, 180, 250, 0.12);
    border-radius: 10px;
    gridline-color: rgba(137, 180, 250, 0.05);
    font-size: 9pt;
    selection-background-color: rgba(137, 180, 250, 0.18);
    selection-color: #dde6ff;
    alternate-background-color: rgba(30, 32, 48, 0.5);
}

QTableView::item, QTableWidget::item {
    padding: 7px 12px;
    color: #c5cef0;
    border-bottom: 1px solid rgba(137, 180, 250, 0.04);
}

QTableView::item:selected, QTableWidget::item:selected {
    background: rgba(137, 180, 250, 0.18);
    color: #dde6ff;
}

QTableView::item:hover, QTableWidget::item:hover {
    background: rgba(137, 180, 250, 0.08);
}

QHeaderView {
    background: transparent;
}

QHeaderView::section {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 rgba(40, 42, 62, 1.0), stop:1 rgba(32, 34, 52, 1.0));
    color: #5a5c78;
    padding: 9px 12px;
    border: none;
    border-right: 1px solid rgba(137, 180, 250, 0.06);
    border-bottom: 1px solid rgba(137, 180, 250, 0.1);
    font-weight: 700;
    font-size: 8pt;
    text-transform: uppercase;
    letter-spacing: 0.07em;
}

QHeaderView::section:first {
    border-top-left-radius: 10px;
}

QHeaderView::section:last {
    border-top-right-radius: 10px;
    border-right: none;
}

QHeaderView::section:hover {
    background: rgba(50, 52, 72, 1.0);
    color: #8a8cb0;
}

/* ════════════════════════════════════════════════════════════════════
   LISTWIDGET — marqueurs, clusters
   ════════════════════════════════════════════════════════════════════ */

QListWidget {
    background: rgba(14, 15, 24, 0.95);
    border: 1px solid rgba(137, 180, 250, 0.12);
    border-radius: 10px;
    padding: 4px;
    font-size: 10pt;
    outline: none;
}

QListWidget::item {
    padding: 7px 12px;
    border-radius: 6px;
    margin: 1px 3px;
    color: #8a8cb0;
}

QListWidget::item:selected {
    background: rgba(137, 180, 250, 0.18);
    color: #dde6ff;
    border: 1px solid rgba(137, 180, 250, 0.28);
}

QListWidget::item:hover:!selected {
    background: rgba(137, 180, 250, 0.07);
    color: #b0b8d8;
}

/* ════════════════════════════════════════════════════════════════════
   SPLITTER
   ════════════════════════════════════════════════════════════════════ */

QSplitter::handle {
    background: rgba(137, 180, 250, 0.08);
    width: 1px;
}

QSplitter::handle:hover {
    background: rgba(137, 180, 250, 0.35);
}

/* ════════════════════════════════════════════════════════════════════
   TOOLTIPS
   ════════════════════════════════════════════════════════════════════ */

QToolTip {
    background: #252637;
    color: #c5cef0;
    border: 1px solid rgba(137, 180, 250, 0.22);
    border-radius: 7px;
    padding: 7px 12px;
    font-size: 9pt;
    font-family: 'Segoe UI', 'Inter', Arial, sans-serif;
}

/* ════════════════════════════════════════════════════════════════════
   FRAMELINE séparateur
   ════════════════════════════════════════════════════════════════════ */

QFrame[frameShape="4"],
QFrame[frameShape="5"] {
    color: rgba(137, 180, 250, 0.1);
    border: none;
    max-height: 1px;
    min-height: 1px;
}

/* ════════════════════════════════════════════════════════════════════
   MESSAGEBOX
   ════════════════════════════════════════════════════════════════════ */

QMessageBox {
    background: #1a1a2e;
    color: #cdd6f4;
}

QMessageBox QLabel {
    color: #cdd6f4;
    font-size: 10pt;
}

QMessageBox QPushButton {
    min-width: 80px;
    padding: 8px 16px;
}
"""
