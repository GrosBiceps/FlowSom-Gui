# -*- coding: utf-8 -*-
"""
styles.py — PRISMA v2.0 Brand Identity System.

Flat design, technical, scientific. Sharp corners, no background gradients.

Palette PRISMA:
    Void        #04070D    Deep      #080D18    Surface   #0C1220
    Raised      #101825    Elevated  #141E2E    Paper     #EEF2F7

    Channels (cytometric spectrum):
      APC       #FF3D6E   (MRD, danger)
      PE        #FF9B3D   (granulocytes, warm)
      PerCP     #FFE032   (caution, labels)
      FITC      #39FF8A   (accent, success)
      V500      #5BAAFF   (T-cells, info)
      V450      #7B52FF   (brand primary, lymphocytes)
      SSC       #7EC8E3   (granularity, steel)
      FSC       #E8F4FF   (size, light)

Typography:
  - UI / display : Outfit → Segoe UI → Inter
  - Data / code  : JetBrains Mono → Cascadia Code → Consolas
  - Radius       : 0px (panels), 2px max (small elements)
"""

# ---------------------------------------------------------------------------
# PRISMA v2.0 — Design Tokens
# ---------------------------------------------------------------------------
COLORS = {
    # Core surfaces
    "void": "#04070D",
    "deep": "#080D18",
    "surface": "#0C1220",
    "raised": "#101825",
    "elevated": "#141E2E",
    "paper": "#EEF2F7",
    "paper_dim": "rgba(238,242,247,0.55)",
    # Borders
    "border": "rgba(255,255,255,0.055)",
    "border_mid": "rgba(255,255,255,0.10)",
    # Cytometric channels
    "ch_apc": "#FF3D6E",
    "ch_pe": "#FF9B3D",
    "ch_percp": "#FFE032",
    "ch_fitc": "#39FF8A",
    "ch_v500": "#5BAAFF",
    "ch_v450": "#7B52FF",
    "ch_ssc": "#7EC8E3",
    "ch_fsc": "#E8F4FF",
    # Brand primaries (semantic aliases)
    "brand": "#7B52FF",
    "accent": "#39FF8A",
    "warm": "#FF9B3D",
    "danger": "#FF3D6E",
    "info": "#5BAAFF",
    "warn": "#FFE032",
}

# Stepper sémantiques — PRISMA channels
STEP_COLORS = {
    "done": COLORS["accent"],  # FITC
    "active": COLORS["brand"],  # V450
    "pending": "#2A3342",  # dim neutral
    "error": COLORS["danger"],  # APC
}

# ---------------------------------------------------------------------------
# QSS — PRISMA v2.0
# ---------------------------------------------------------------------------
STYLESHEET = """
/* ════════════════════════════════════════════════════════════════════
   BASE — Void foundation
   ════════════════════════════════════════════════════════════════════ */

QMainWindow {
    background: #04070D;
}

QWidget {
    background-color: transparent;
    color: #EEF2F7;
    font-family: 'Segoe UI', Arial, sans-serif;
    font-size: 10pt;
    font-weight: 400;
}

QWidget#centralWidget,
QWidget#wizardShell {
    background: #04070D;
}

/* Sidebar — Deep */
QWidget#stepSidebar,
QWidget#sidebarHeader {
    background: #080D18;
    border-right: 1px solid rgba(255,255,255,0.055);
    border-bottom: none;
}

/* Main content — Surface */
QWidget#stepContent,
QWidget#rightPanel {
    background: #0C1220;
    border-radius: 0px;
}

QWidget#leftPanel {
    background: #080D18;
    border-right: 1px solid rgba(255,255,255,0.055);
}

/* ════════════════════════════════════════════════════════════════════
   SIDEBAR HEADER — branding PRISMA
   ════════════════════════════════════════════════════════════════════ */

QLabel#brandTitle {
    color: #EEF2F7;
    font-family: 'Segoe UI', sans-serif;
    font-size: 18pt;
    font-weight: 900;
    letter-spacing: -0.04em;
    background: transparent;
}

QLabel#brandSubtitle {
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 8pt;
    font-weight: 400;
    letter-spacing: 0.22em;
    background: transparent;
}

QLabel#sidebarSeparator,
QFrame#sidebarSeparator {
    background: rgba(255,255,255,0.055);
    color: rgba(255,255,255,0.055);
    border: none;
    max-height: 1px;
    min-height: 1px;
}

QLabel#versionLabel {
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 8pt;
    letter-spacing: 0.08em;
    background: transparent;
    padding: 10px 18px;
}

/* ════════════════════════════════════════════════════════════════════
   STEPPER — navigation latérale (flat, border-left)
   ════════════════════════════════════════════════════════════════════ */

QLabel#stepLabel {
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 8pt;
    font-weight: 500;
    letter-spacing: 0.22em;
    text-transform: uppercase;
    padding: 3px 0px;
}

QLabel#stepTitle {
    font-size: 10pt;
    font-weight: 700;
    color: #EEF2F7;
}

QLabel#stepDone    { color: #39FF8A; }
QLabel#stepActive  { color: #7B52FF; }
QLabel#stepPending { color: #2A3342; }

QPushButton#stepBtn {
    background: transparent;
    border: none;
    border-left: 2px solid transparent;
    border-radius: 0px;
    padding: 12px 18px;
    text-align: left;
    font-family: 'Segoe UI', sans-serif;
    font-size: 10pt;
    font-weight: 400;
    color: #EEF2F7;
}

QPushButton#stepBtn:hover {
    background: rgba(123,82,255,0.06);
    color: #EEF2F7;
    border-left-color: rgba(123,82,255,0.35);
}

QPushButton#stepBtn:focus {
    background: rgba(123,82,255,0.10);
    color: #EEF2F7;
    border-left-color: rgba(123,82,255,0.70);
    outline: none;
}

QPushButton#stepBtnActive {
    background: rgba(123,82,255,0.12);
    border: none;
    border-left: 2px solid #7B52FF;
    border-radius: 0px;
    padding: 12px 18px;
    text-align: left;
    font-family: 'Segoe UI', sans-serif;
    font-size: 10pt;
    font-weight: 700;
    color: #EEF2F7;
}

QPushButton#stepBtnDone {
    background: transparent;
    border: none;
    border-left: 2px solid rgba(57,255,138,0.55);
    border-radius: 0px;
    padding: 12px 18px;
    text-align: left;
    font-family: 'Segoe UI', sans-serif;
    font-size: 10pt;
    font-weight: 500;
    color: rgba(57,255,138,0.85);
}

QPushButton#stepBtnDone:hover {
    background: rgba(57,255,138,0.06);
    color: #39FF8A;
}

QPushButton#stepBtnDone:focus {
    background: rgba(57,255,138,0.12);
    color: #39FF8A;
    border-left-color: #39FF8A;
    outline: none;
}

/* ════════════════════════════════════════════════════════════════════
   TITLE BARS — entêtes plates PRISMA
   ════════════════════════════════════════════════════════════════════ */

QWidget#titleBar,
QWidget#actionBar,
QWidget#navBar {
    background: #080D18;
    border-bottom: 1px solid rgba(255,255,255,0.055);
}

QWidget#navBar {
    border-bottom: none;
    border-top: 1px solid rgba(255,255,255,0.055);
}

/* ════════════════════════════════════════════════════════════════════
   GROUPBOX — cartes flat (Raised)
   ════════════════════════════════════════════════════════════════════ */

QGroupBox {
    font-family: 'Segoe UI', sans-serif;
    font-weight: 700;
    font-size: 10pt;
    letter-spacing: 0.02em;
    border: 1px solid rgba(255,255,255,0.055);
    border-top: 2px solid rgba(123,82,255,0.28);
    border-radius: 0px;
    margin-top: 18px;
    padding: 18px 14px 14px 14px;
    background: #101825;
}

QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 14px;
    top: 2px;
    padding: 2px 10px;
    color: #EEF2F7;
    background: #101825;
    border: none;
    border-radius: 0px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    font-weight: 500;
    letter-spacing: 0.18em;
    text-transform: uppercase;
}

/* ════════════════════════════════════════════════════════════════════
   DRAG & DROP — Import
   ════════════════════════════════════════════════════════════════════ */

QLabel#dropZone {
    background: rgba(123,82,255,0.03);
    border: 1px dashed rgba(123,82,255,0.30);
    border-radius: 0px;
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 10pt;
    padding: 32px 20px;
    qproperty-alignment: AlignCenter;
}

QLabel#dropZone[dragOver="true"] {
    background: rgba(123,82,255,0.10);
    border: 1px dashed #7B52FF;
    color: #EEF2F7;
}

QLabel#dropZoneOk {
    background: rgba(57,255,138,0.05);
    border: 1px solid rgba(57,255,138,0.38);
    border-radius: 0px;
    color: #39FF8A;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 10pt;
    padding: 32px 20px;
    qproperty-alignment: AlignCenter;
}

/* ════════════════════════════════════════════════════════════════════
   BOUTONS — PRISMA flat (btn-fill / outline / ghost / danger)
   ════════════════════════════════════════════════════════════════════ */

/* Base = btn-ghost (secondaire) */
QPushButton {
    background: rgba(255,255,255,0.04);
    border: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
    padding: 9px 18px;
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    font-weight: 500;
    letter-spacing: 0.10em;
    text-transform: uppercase;
    min-height: 14px;
}

QPushButton:hover {
    background: rgba(255,255,255,0.08);
    border-color: rgba(255,255,255,0.10);
    color: #EEF2F7;
}

QPushButton:focus {
    border-color: rgba(123,82,255,0.70);
    background: rgba(123,82,255,0.12);
    color: #EEF2F7;
    outline: none;
}

QPushButton:pressed {
    background: rgba(255,255,255,0.03);
}

QPushButton:disabled {
    background: rgba(255,255,255,0.02);
    color: #EEF2F7;
    border-color: rgba(255,255,255,0.03);
}

/* PRIMAIRE — btn-fill (brand violet) */
QPushButton#primaryBtn {
    background: #7B52FF;
    border: 1px solid #7B52FF;
    border-radius: 0px;
    color: #FFFFFF;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 700;
    font-size: 9pt;
    letter-spacing: 0.14em;
    text-transform: uppercase;
    padding: 11px 26px;
}

QPushButton#primaryBtn:hover {
    background: #8F68FF;
    border-color: #8F68FF;
    color: #FFFFFF;
}

QPushButton#primaryBtn:focus {
    border-color: #A98EFF;
    background: #8F68FF;
}

QPushButton#primaryBtn:pressed {
    background: #6A42E6;
    border-color: #6A42E6;
}

QPushButton#primaryBtn:disabled {
    background: rgba(123,82,255,0.35);
    border-color: rgba(123,82,255,0.35);
    color: rgba(255,255,255,0.55);
}

/* SUCCÈS — accent vert FITC */
QPushButton#successBtn {
    background: transparent;
    border: 1px solid rgba(57,255,138,0.55);
    border-radius: 0px;
    color: #39FF8A;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 700;
    font-size: 9pt;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    padding: 9px 18px;
}

QPushButton#successBtn:hover {
    background: rgba(57,255,138,0.10);
    border-color: #39FF8A;
    color: #39FF8A;
}

QPushButton#successBtn:focus {
    background: rgba(57,255,138,0.16);
    border-color: #39FF8A;
}

QPushButton#successBtn:pressed {
    background: rgba(57,255,138,0.18);
}

QPushButton#successBtn:disabled {
    background: transparent;
    color: rgba(57,255,138,0.25);
    border-color: rgba(57,255,138,0.18);
}

/* DANGER — APC (rouge) */
QPushButton#dangerBtn {
    background: rgba(255,61,110,0.12);
    border: 1px solid rgba(255,61,110,0.35);
    border-radius: 0px;
    color: #FF3D6E;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 700;
    font-size: 9pt;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    padding: 9px 18px;
}

QPushButton#dangerBtn:hover {
    background: rgba(255,61,110,0.22);
    border-color: #FF3D6E;
    color: #FF3D6E;
}

QPushButton#dangerBtn:focus {
    background: rgba(255,61,110,0.26);
    border-color: #FF3D6E;
}

QPushButton#dangerBtn:pressed {
    background: rgba(255,61,110,0.30);
}

QPushButton#dangerBtn:disabled {
    background: rgba(255,61,110,0.05);
    color: rgba(255,61,110,0.25);
    border-color: rgba(255,61,110,0.12);
}

/* EXPORT — V500 (info) — traité en outline info */
QPushButton#exportBtn {
    background: transparent;
    border: 1px solid rgba(91,170,255,0.32);
    border-radius: 0px;
    color: #5BAAFF;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 500;
    font-size: 9pt;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    padding: 8px 16px;
}

QPushButton#exportBtn:hover {
    background: rgba(91,170,255,0.10);
    border-color: #5BAAFF;
    color: #5BAAFF;
}

QPushButton#exportBtn:focus {
    background: rgba(91,170,255,0.16);
    border-color: #5BAAFF;
}

QPushButton#exportBtn:pressed {
    background: rgba(91,170,255,0.18);
}

QPushButton#exportBtn:disabled {
    background: transparent;
    color: rgba(91,170,255,0.22);
    border-color: rgba(91,170,255,0.12);
}

/* GHOST — outline fin (btn-outline HTML) */
QPushButton#ghostBtn {
    background: transparent;
    border: 1px solid rgba(255,255,255,0.10);
    border-radius: 0px;
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 500;
    font-size: 9pt;
    letter-spacing: 0.06em;
    padding: 8px 14px;
}

QPushButton#ghostBtn:hover {
    background: transparent;
    border-color: #39FF8A;
    color: #39FF8A;
}

QPushButton#ghostBtn:focus {
    background: rgba(57,255,138,0.08);
    border-color: #39FF8A;
    color: #39FF8A;
}

QPushButton#ghostBtn:pressed {
    background: rgba(57,255,138,0.08);
    border-color: #39FF8A;
}

/* KEEP — identique à successBtn */
QPushButton#keepBtn {
    background: rgba(57,255,138,0.08);
    border: 1px solid #39FF8A;
    border-radius: 0px;
    color: #39FF8A;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 700;
    font-size: 9pt;
    letter-spacing: 0.14em;
    text-transform: uppercase;
    padding: 10px 26px;
    min-height: 28px;
}

QPushButton#keepBtn:hover {
    background: rgba(57,255,138,0.18);
    border-color: #39FF8A;
}

QPushButton#keepBtn:focus {
    background: rgba(57,255,138,0.22);
    border-color: #39FF8A;
}

QPushButton#keepBtn:pressed {
    background: rgba(57,255,138,0.28);
}

/* DISCARD — identique à dangerBtn */
QPushButton#discardBtn {
    background: rgba(255,61,110,0.10);
    border: 1px solid #FF3D6E;
    border-radius: 0px;
    color: #FF3D6E;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 700;
    font-size: 9pt;
    letter-spacing: 0.14em;
    text-transform: uppercase;
    padding: 10px 26px;
    min-height: 28px;
}

QPushButton#discardBtn:hover {
    background: rgba(255,61,110,0.22);
    border-color: #FF3D6E;
}

QPushButton#discardBtn:focus {
    background: rgba(255,61,110,0.28);
    border-color: #FF3D6E;
}

QPushButton#discardBtn:pressed {
    background: rgba(255,61,110,0.32);
}

/* ICON-ONLY */
QPushButton#iconBtn {
    background: rgba(255,255,255,0.03);
    border: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
    padding: 6px;
    min-width: 32px;
    max-width: 36px;
    min-height: 32px;
    max-height: 36px;
}

QPushButton#iconBtn:hover {
    background: rgba(123,82,255,0.12);
    border-color: rgba(123,82,255,0.35);
}

QPushButton#iconBtn:focus {
    background: rgba(123,82,255,0.18);
    border-color: rgba(123,82,255,0.60);
}

/* ════════════════════════════════════════════════════════════════════
   LABELS
   ════════════════════════════════════════════════════════════════════ */

QLabel {
    color: #EEF2F7;
    background: transparent;
    font-family: 'Segoe UI', sans-serif;
}

QLabel#titleLabel {
    font-family: 'Segoe UI', sans-serif;
    font-size: 22pt;
    font-weight: 900;
    color: #EEF2F7;
    letter-spacing: -0.04em;
    background: transparent;
}

QLabel#subtitleLabel {
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    color: #EEF2F7;
    font-weight: 400;
    letter-spacing: 0.04em;
    background: transparent;
}

QLabel#sectionLabel {
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    font-weight: 500;
    color: #39FF8A;
    letter-spacing: 0.22em;
    text-transform: uppercase;
    background: transparent;
}

QLabel#cardLabel {
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 8pt;
    color: #EEF2F7;
    font-weight: 500;
    text-transform: uppercase;
    letter-spacing: 0.18em;
    background: transparent;
}

QLabel#fileLabel {
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    color: #EEF2F7;
    padding: 7px 12px;
    background: rgba(255,255,255,0.03);
    border: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
}

QLabel#hintLabel {
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    letter-spacing: 0.04em;
    background: transparent;
}

QLabel#summaryLabel {
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    padding: 4px 2px;
    background: transparent;
}

QLabel#summaryLabelActive {
    color: #7B52FF;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    font-weight: 700;
    padding: 4px 2px;
    background: transparent;
}

/* States — success / error / warning */
QLabel#statusSuccess {
    color: #39FF8A;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 500;
}

QLabel#statusError {
    color: #FF3D6E;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 500;
}

QLabel#statusWarning {
    color: #FFE032;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 500;
}

/* Badge pill (ch-badge PRISMA) */
QLabel#badge {
    background: rgba(123,82,255,0.18);
    color: #7B52FF;
    border-radius: 0px;
    border: 1px solid rgba(123,82,255,0.32);
    padding: 3px 10px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 8pt;
    font-weight: 500;
    letter-spacing: 0.10em;
    text-transform: uppercase;
}

QLabel#badgeGreen {
    background: rgba(57,255,138,0.12);
    color: #39FF8A;
    border-radius: 0px;
    border: 1px solid rgba(57,255,138,0.28);
    padding: 3px 10px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 8pt;
    font-weight: 500;
    letter-spacing: 0.10em;
    text-transform: uppercase;
}

QLabel#badgeRed {
    background: rgba(255,61,110,0.15);
    color: #FF3D6E;
    border-radius: 0px;
    border: 1px solid rgba(255,61,110,0.30);
    padding: 3px 10px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 8pt;
    font-weight: 500;
    letter-spacing: 0.10em;
    text-transform: uppercase;
}

/* ════════════════════════════════════════════════════════════════════
   ALERTS — bandeaux sémantiques (HTML .alert-*)
   ════════════════════════════════════════════════════════════════════ */

QLabel#alertSuccess, QWidget#alertSuccess {
    border: 1px solid rgba(57,255,138,0.30);
    color: rgba(57,255,138,0.85);
    background: rgba(57,255,138,0.05);
    padding: 10px 16px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 10pt;
    border-radius: 0px;
}

QLabel#alertWarn, QWidget#alertWarn {
    border: 1px solid rgba(255,224,50,0.30);
    color: rgba(255,224,50,0.85);
    background: rgba(255,224,50,0.05);
    padding: 10px 16px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 10pt;
    border-radius: 0px;
}

QLabel#alertDanger, QWidget#alertDanger {
    border: 1px solid rgba(255,61,110,0.30);
    color: rgba(255,61,110,0.85);
    background: rgba(255,61,110,0.05);
    padding: 10px 16px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 10pt;
    border-radius: 0px;
}

QLabel#alertInfo, QWidget#alertInfo {
    border: 1px solid rgba(91,170,255,0.30);
    color: rgba(91,170,255,0.85);
    background: rgba(91,170,255,0.05);
    padding: 10px 16px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 10pt;
    border-radius: 0px;
}

/* Bandeau d'avertissement clinique (outil de recherche) */
QWidget#clinicalWarningBanner {
    background: rgba(255,61,110,0.06);
    border-bottom: 1px solid rgba(255,61,110,0.30);
}

QLabel#clinicalWarningText {
    background: transparent;
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
}

/* ════════════════════════════════════════════════════════════════════
   VALIDATION
   ════════════════════════════════════════════════════════════════════ */

QSpinBox[invalid="true"],
QDoubleSpinBox[invalid="true"],
QLineEdit[invalid="true"] {
    border: 1px solid #FF3D6E;
    background: rgba(255,61,110,0.07);
}

/* ════════════════════════════════════════════════════════════════════
   INPUTS — flat, monospace
   ════════════════════════════════════════════════════════════════════ */

QSpinBox, QDoubleSpinBox, QComboBox, QLineEdit {
    background: rgba(255,255,255,0.03);
    border: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
    padding: 8px 10px;
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 10pt;
    min-height: 14px;
    selection-background-color: rgba(123,82,255,0.35);
    selection-color: #EEF2F7;
}

QSpinBox:hover, QDoubleSpinBox:hover, QComboBox:hover, QLineEdit:hover {
    border-color: rgba(255,255,255,0.10);
    background: rgba(255,255,255,0.05);
}

QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus, QLineEdit:focus {
    border: 1px solid rgba(123,82,255,0.50);
    background: rgba(255,255,255,0.04);
    outline: none;
}

QSpinBox:disabled, QDoubleSpinBox:disabled, QComboBox:disabled, QLineEdit:disabled {
    color: #EEF2F7;
    border-color: rgba(255,255,255,0.03);
    background: rgba(255,255,255,0.02);
}

QLineEdit[state="success"] {
    border-color: rgba(57,255,138,0.40);
}

QLineEdit[state="error"] {
    border-color: rgba(255,61,110,0.40);
}

QSpinBox::up-button, QSpinBox::down-button,
QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {
    background: rgba(123,82,255,0.08);
    border: none;
    border-radius: 0px;
    width: 16px;
    margin: 2px;
}

QSpinBox::up-button:hover, QSpinBox::down-button:hover,
QDoubleSpinBox::up-button:hover, QDoubleSpinBox::down-button:hover {
    background: rgba(123,82,255,0.22);
}

QComboBox::drop-down {
    border: none;
    width: 24px;
    background: rgba(123,82,255,0.08);
    border-radius: 0px;
}

QComboBoxPrivateContainer {
    background: #101825;
    border: 1px solid rgba(255,255,255,0.10);
    margin: 0px;
    padding: 0px;
}

QComboBoxPrivateScroller {
    background: #101825;
    border: none;
    margin: 0px;
    padding: 0px;
    min-height: 0px;
}

QComboBoxPrivateContainer QWidget {
    background: #101825;
    border: none;
}

QComboBox QAbstractItemView {
    background: #101825;
    border: 1px solid #101825;
    selection-background-color: rgba(123,82,255,0.24);
    selection-color: #EEF2F7;
    outline: none;
    padding: 2px;
}

QComboBox QAbstractItemView::item {
    background: #101825;
    color: #EEF2F7;
    padding: 7px 12px;
    border: none;
    min-height: 20px;
    border-radius: 0px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
}

QComboBox QAbstractItemView::item:hover,
QComboBox QAbstractItemView::item:focus {
    background: rgba(123,82,255,0.14);
}

QComboBox QAbstractItemView::item:selected {
    background: rgba(123,82,255,0.28);
    color: #EEF2F7;
}

QComboBox QScrollBar:vertical {
    background: #101825;
    width: 6px;
    border: none;
}

QComboBox QScrollBar::handle:vertical {
    background: rgba(123,82,255,0.35);
    border-radius: 0px;
    min-height: 20px;
}

QComboBox QFrame {
    background: #101825;
    border: none;
}

/* ════════════════════════════════════════════════════════════════════
   CHECKBOX / RADIOBUTTON
   ════════════════════════════════════════════════════════════════════ */

QCheckBox {
    spacing: 10px;
    color: #EEF2F7;
    font-family: 'Segoe UI', sans-serif;
    font-size: 10pt;
}

QCheckBox::indicator {
    width: 16px;
    height: 16px;
    border-radius: 0px;
    border: 1px solid rgba(255,255,255,0.18);
    background: rgba(255,255,255,0.03);
}

QCheckBox::indicator:hover {
    border-color: rgba(123,82,255,0.55);
    background: rgba(123,82,255,0.08);
}

QCheckBox::indicator:focus {
    border-color: rgba(123,82,255,0.70);
    background: rgba(123,82,255,0.14);
}

QCheckBox::indicator:checked {
    background: #7B52FF;
    border-color: #7B52FF;
    image: none;
}

QCheckBox::indicator:checked:hover {
    background: #8F68FF;
    border-color: #8F68FF;
}

QCheckBox:disabled {
    color: #EEF2F7;
}

QCheckBox::indicator:disabled {
    border-color: rgba(255,255,255,0.04);
    background: rgba(255,255,255,0.02);
}

/* ════════════════════════════════════════════════════════════════════
   PROGRESSBAR — 2-4px fine, Brand → Accent
   ════════════════════════════════════════════════════════════════════ */

QProgressBar {
    border: none;
    border-radius: 0px;
    background: rgba(255,255,255,0.07);
    text-align: center;
    color: transparent;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 500;
    font-size: 8pt;
    min-height: 3px;
    max-height: 3px;
}

QProgressBar::chunk {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 #7B52FF, stop:1 #39FF8A);
    border-radius: 0px;
}

/* Main pipeline progress — 4px fin */
QProgressBar#pipelineProgress {
    border: none;
    border-radius: 0px;
    background: rgba(255,255,255,0.07);
    text-align: center;
    color: transparent;
    min-height: 4px;
    max-height: 4px;
}

QProgressBar#pipelineProgress::chunk {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 #7B52FF, stop:1 #39FF8A);
    border-radius: 0px;
}

/* Step label sous la progress bar (HTML .alert-*) */
QLabel#pipelineStepLabel {
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 10pt;
    font-weight: 400;
    letter-spacing: 0.04em;
    padding: 6px 0;
}

QLabel#pipelineStepLabel[running="true"] {
    color: #7B52FF;
}

QLabel#pipelineStepLabel[done="true"] {
    color: #39FF8A;
}

QLabel#pipelineStepLabel[error="true"] {
    color: #FF3D6E;
}

/* ════════════════════════════════════════════════════════════════════
   LOG CONSOLE — terminal PRISMA
   ════════════════════════════════════════════════════════════════════ */

QTextEdit#logConsole,
QPlainTextEdit#logConsole {
    background: rgba(0,0,0,0.30);
    border: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
    color: #EEF2F7;
    padding: 16px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    font-weight: 300;
    selection-background-color: rgba(123,82,255,0.28);
}

QTextEdit {
    background: rgba(0,0,0,0.25);
    border: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
    color: #EEF2F7;
    padding: 12px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    selection-background-color: rgba(123,82,255,0.28);
}

/* ════════════════════════════════════════════════════════════════════
   STATUSBAR
   ════════════════════════════════════════════════════════════════════ */

QStatusBar {
    background: #04070D;
    color: #EEF2F7;
    border-top: 1px solid rgba(255,255,255,0.055);
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    letter-spacing: 0.08em;
    padding: 4px 12px;
}

/* ════════════════════════════════════════════════════════════════════
   TABWIDGET
   ════════════════════════════════════════════════════════════════════ */

QTabWidget::pane {
    border: 1px solid rgba(255,255,255,0.055);
    border-top: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
    background: #0C1220;
    padding-top: 0px;
    top: 0px;
}

QTabWidget::tab-bar {
    left: 0px;
    alignment: left;
}

QTabBar {
    qproperty-drawBase: 0;
    background: #080D18;
}

QTabBar::tab {
    background: transparent;
    border: none;
    border-bottom: 2px solid transparent;
    padding: 8px 14px 8px 10px;
    margin-right: 2px;
    color: #EEF2F7;
    font-family: 'Segoe UI', sans-serif;
    font-weight: 400;
    font-size: 9pt;
    min-width: 80px;
    text-align: left;
}

QTabBar::tab:selected {
    background: transparent;
    color: #39FF8A;
    font-weight: 500;
    border-bottom: 2px solid #39FF8A;
}

QTabBar::tab:hover:!selected {
    color: #EEF2F7;
    border-bottom: 2px solid rgba(57,255,138,0.25);
}

/* ════════════════════════════════════════════════════════════════════
   SCROLLBAR — 6px PRISMA
   ════════════════════════════════════════════════════════════════════ */

QScrollArea {
    border: none;
    background: transparent;
}

QScrollBar:vertical {
    background: #04070D;
    width: 6px;
    border-radius: 0px;
    margin: 0px;
}

QScrollBar::handle:vertical {
    background: rgba(123,82,255,0.35);
    border-radius: 0px;
    min-height: 28px;
}

QScrollBar::handle:vertical:hover {
    background: rgba(123,82,255,0.55);
}

QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0;
}

QScrollBar:horizontal {
    background: #04070D;
    height: 6px;
    border-radius: 0px;
    margin: 0px;
}

QScrollBar::handle:horizontal {
    background: rgba(123,82,255,0.35);
    border-radius: 0px;
    min-width: 28px;
}

QScrollBar::handle:horizontal:hover {
    background: rgba(123,82,255,0.55);
}

QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
    width: 0;
}

/* ════════════════════════════════════════════════════════════════════
   TABLE
   ════════════════════════════════════════════════════════════════════ */

QTableView, QTableWidget {
    background: #101825;
    border: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
    gridline-color: rgba(255,255,255,0.03);
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    selection-background-color: rgba(123,82,255,0.18);
    selection-color: #EEF2F7;
    alternate-background-color: #141E2E;
}

QTableView::item, QTableWidget::item {
    padding: 8px 12px;
    color: #EEF2F7;
    border-bottom: 1px solid rgba(255,255,255,0.02);
}

QTableView::item:selected, QTableWidget::item:selected {
    background: rgba(123,82,255,0.18);
    color: #EEF2F7;
}

QTableView::item:hover, QTableWidget::item:hover {
    background: rgba(123,82,255,0.06);
}

QHeaderView {
    background: transparent;
}

QHeaderView::section {
    background: #080D18;
    color: #EEF2F7;
    padding: 9px 12px;
    border: none;
    border-right: 1px solid rgba(255,255,255,0.04);
    border-bottom: 1px solid rgba(255,255,255,0.055);
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 500;
    font-size: 8pt;
    text-transform: uppercase;
    letter-spacing: 0.18em;
}

QHeaderView::section:first,
QHeaderView::section:last {
    border-radius: 0px;
}

QHeaderView::section:last {
    border-right: none;
}

QHeaderView::section:hover {
    background: #101825;
    color: #EEF2F7;
}

/* FCS Preview table */
QTableWidget#fcsPreviewTable {
    background: #101825;
    border: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
    gridline-color: rgba(255,255,255,0.03);
    alternate-background-color: #141E2E;
}

QTableWidget#fcsPreviewTable::item {
    padding: 8px 14px;
    color: #EEF2F7;
    border-bottom: 1px solid rgba(255,255,255,0.02);
}

QTableWidget#fcsPreviewTable::item[status="ok"] { color: #39FF8A; }
QTableWidget#fcsPreviewTable::item[status="warn"] { color: #FFE032; }

QLabel#statusBadgeOk {
    background: rgba(57,255,138,0.12);
    color: #39FF8A;
    border: 1px solid rgba(57,255,138,0.30);
    border-radius: 0px;
    padding: 2px 8px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 8pt;
    font-weight: 500;
    letter-spacing: 0.10em;
    text-transform: uppercase;
}

QLabel#statusBadgePatho {
    background: rgba(255,61,110,0.14);
    color: #FF3D6E;
    border: 1px solid rgba(255,61,110,0.30);
    border-radius: 0px;
    padding: 2px 8px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 8pt;
    font-weight: 500;
    letter-spacing: 0.10em;
    text-transform: uppercase;
}

QLabel#statusBadgeWarn {
    background: rgba(255,224,50,0.12);
    color: #FFE032;
    border: 1px solid rgba(255,224,50,0.30);
    border-radius: 0px;
    padding: 2px 8px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 8pt;
    font-weight: 500;
    letter-spacing: 0.10em;
    text-transform: uppercase;
}

/* ════════════════════════════════════════════════════════════════════
   LISTWIDGET
   ════════════════════════════════════════════════════════════════════ */

QListWidget {
    background: #101825;
    border: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
    padding: 4px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    outline: none;
}

QListWidget::item {
    padding: 7px 12px;
    border-radius: 0px;
    margin: 1px 2px;
    color: #EEF2F7;
}

QListWidget::item:selected {
    background: rgba(123,82,255,0.18);
    color: #EEF2F7;
    border: 1px solid rgba(123,82,255,0.35);
}

QListWidget::item:hover:!selected {
    background: rgba(123,82,255,0.06);
    color: #EEF2F7;
}

/* ════════════════════════════════════════════════════════════════════
   SPLITTER
   ════════════════════════════════════════════════════════════════════ */

QSplitter::handle {
    background: rgba(255,255,255,0.055);
    width: 1px;
}

QSplitter::handle:hover {
    background: rgba(123,82,255,0.45);
}

/* ════════════════════════════════════════════════════════════════════
   TOOLTIPS
   ════════════════════════════════════════════════════════════════════ */

QToolTip {
    background: #141E2E;
    color: #EEF2F7;
    border: 1px solid rgba(255,255,255,0.10);
    border-radius: 0px;
    padding: 7px 12px;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
}

/* ════════════════════════════════════════════════════════════════════
   FRAMES / SEPARATORS
   ════════════════════════════════════════════════════════════════════ */

QFrame[frameShape="4"],
QFrame[frameShape="5"] {
    color: rgba(255,255,255,0.055);
    background: rgba(255,255,255,0.055);
    border: none;
    max-height: 1px;
    min-height: 1px;
}

QFrame#settingsCard {
    background: #101825;
    border: 1px solid rgba(255,255,255,0.055);
    border-top: 2px solid rgba(123,82,255,0.28);
    border-radius: 0px;
}

/* ════════════════════════════════════════════════════════════════════
   WELCOME SCREEN (ÉTAPE 0)
   ════════════════════════════════════════════════════════════════════ */

QWidget#welcomePage {
    background: qradialgradient(cx:0.70, cy:0.22, radius:1.0,
        fx:0.70, fy:0.22,
        stop:0 rgba(123,82,255,0.12),
        stop:0.30 rgba(91,170,255,0.06),
        stop:0.55 rgba(57,255,138,0.03),
        stop:1 #04070D);
}

QWidget#welcomeCard {
    background: #0C1220;
    border: 1px solid rgba(255,255,255,0.055);
    border-top: 1px solid rgba(123,82,255,0.42);
    border-radius: 0px;
}

QFrame#welcomeSpectrum {
    border: none;
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 rgba(255,61,110,0.00),
        stop:0.08 rgba(255,61,110,0.70),
        stop:0.22 rgba(255,155,61,0.70),
        stop:0.37 rgba(255,224,50,0.55),
        stop:0.50 rgba(57,255,138,0.72),
        stop:0.66 rgba(91,170,255,0.72),
        stop:0.82 rgba(123,82,255,0.72),
        stop:1 rgba(232,244,255,0.00));
}

QLabel#welcomeEyebrow {
    color: #39FF8A;
    font-family: 'Segoe UI', sans-serif;
    font-size: 9pt;
    font-weight: 700;
    letter-spacing: 0.18em;
    text-transform: uppercase;
    background: transparent;
    padding-top: 6px;
}

QLabel#welcomeHeroTitle {
    font-family: 'Segoe UI', sans-serif;
    font-size: 31pt;
    font-weight: 900;
    color: #EEF2F7;
    letter-spacing: -0.055em;
    background: transparent;
}

QLabel#welcomeExpansion {
    color: #EEF2F7;
    font-family: 'Segoe UI', sans-serif;
    font-size: 10.2pt;
    font-weight: 300;
    letter-spacing: 0.01em;
    background: transparent;
}

QLabel#welcomeExpansion b {
    color: #EEF2F7;
}

QWidget#welcomeChannelRow {
    background: transparent;
}

QLabel#welcomeChannelBadge {
    border: 1px solid rgba(255,255,255,0.10);
    background: rgba(255,255,255,0.03);
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 7.4pt;
    font-weight: 500;
    letter-spacing: 0.10em;
    text-transform: uppercase;
    padding: 4px 10px;
}

QLabel#welcomeChannelBadge[channel="fitc"] {
    color: #39FF8A;
    border-color: rgba(57,255,138,0.35);
    background: rgba(57,255,138,0.08);
}

QLabel#welcomeChannelBadge[channel="pe"] {
    color: #FF9B3D;
    border-color: rgba(255,155,61,0.35);
    background: rgba(255,155,61,0.08);
}

QLabel#welcomeChannelBadge[channel="apc"] {
    color: #FF3D6E;
    border-color: rgba(255,61,110,0.35);
    background: rgba(255,61,110,0.08);
}

QLabel#welcomeChannelBadge[channel="v450"] {
    color: #7B52FF;
    border-color: rgba(123,82,255,0.35);
    background: rgba(123,82,255,0.08);
}

QWidget#welcomeBody {
    background: #0F1728;
    border: 1px solid rgba(255,255,255,0.05);
}

QWidget#welcomeColLeft,
QWidget#welcomeColRight {
    background: #101825;
    border: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
}

QWidget#welcomeColLeft {
    border-top: 1px solid rgba(57,255,138,0.35);
}

QWidget#welcomeColRight {
    border-top: 1px solid rgba(91,170,255,0.35);
}

QLabel#welcomeColTitle {
    color: #EEF2F7;
    font-family: 'Segoe UI', sans-serif;
    font-size: 9pt;
    font-weight: 700;
    background: transparent;
    letter-spacing: 0.16em;
    text-transform: uppercase;
    padding-bottom: 3px;
}

QLabel#welcomeColItem {
    color: #EEF2F7;
    font-family: 'Segoe UI', sans-serif;
    font-size: 9.6pt;
    font-weight: 400;
    background: transparent;
    letter-spacing: 0.01em;
    padding: 2px 0px;
}

QLabel#welcomeColItem[lane="left"] {
    color: #39FF8A;
}

QLabel#welcomeColItem[lane="right"] {
    color: #5BAAFF;
}

QLabel#welcomeSub {
    color: #EEF2F7;
    font-family: 'Georgia', 'Times New Roman', serif;
    font-style: italic;
    font-size: 12.5pt;
    background: transparent;
}

QLabel#welcomeInfo {
    color: #EEF2F7;
    font-family: 'Segoe UI', sans-serif;
    font-size: 9.4pt;
    font-weight: 300;
    background: transparent;
    letter-spacing: 0.01em;
}

QWidget#welcomeStatsRow {
    background: transparent;
}

QWidget#welcomeStatCard {
    background: #101825;
    border: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
    border-top: 1px solid rgba(255,255,255,0.08);
}

QLabel#welcomeStatValue {
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 13pt;
    font-weight: 700;
    color: #EEF2F7;
    letter-spacing: -0.02em;
    background: transparent;
}

QLabel#welcomeStatValue[tone="accent"] { color: #39FF8A; }
QLabel#welcomeStatValue[tone="brand"] { color: #7B52FF; }
QLabel#welcomeStatValue[tone="warm"] { color: #FF9B3D; }
QLabel#welcomeStatValue[tone="info"] { color: #5BAAFF; }

QLabel#welcomeStatKey {
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 7.8pt;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    background: transparent;
}

QPushButton#welcomeCta {
    background: #7B52FF;
    border: 1px solid #7B52FF;
    border-radius: 0px;
    color: #FFFFFF;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-weight: 700;
    font-size: 9pt;
    letter-spacing: 0.14em;
    text-transform: uppercase;
    padding: 12px 30px;
}

QPushButton#welcomeCta:hover {
    background: #8F68FF;
    border-color: #8F68FF;
    color: #FFFFFF;
}

QPushButton#welcomeCta:pressed {
    background: #6A42E6;
    border-color: #6A42E6;
}

QLabel#welcomeHint {
    color: #EEF2F7;
    font-family: 'Segoe UI', sans-serif;
    font-size: 9pt;
    font-weight: 300;
    letter-spacing: 0.02em;
    background: transparent;
}

/* ════════════════════════════════════════════════════════════════════
   DIALOGS & MESSAGEBOX
   ════════════════════════════════════════════════════════════════════ */

QDialog {
    background: #080D18;
    color: #EEF2F7;
}

QLabel#dialogTitle {
    color: #EEF2F7;
    font-family: 'Segoe UI', sans-serif;
    font-size: 14pt;
    font-weight: 900;
    letter-spacing: -0.03em;
    background: transparent;
}

QLabel#dialogDesc {
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    padding-bottom: 4px;
    background: transparent;
}

QLabel#dialogCount {
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    background: transparent;
    letter-spacing: 0.10em;
}

QMessageBox {
    background: #080D18;
    color: #EEF2F7;
}

QMessageBox QLabel {
    color: #EEF2F7;
    font-family: 'Segoe UI', sans-serif;
    font-size: 10pt;
}

QMessageBox QPushButton {
    min-width: 80px;
    padding: 8px 16px;
}

/* ════════════════════════════════════════════════════════════════════
   MATPLOTLIB TOOLBAR
   ════════════════════════════════════════════════════════════════════ */

QToolBar#matplotlibToolbar {
    background: #0C1220;
    border: 1px solid rgba(255,255,255,0.055);
    border-radius: 0px;
    padding: 3px;
}

/* FCS info label */
QLabel#fcsInfoLabel {
    color: #EEF2F7;
    font-family: 'Consolas', 'Cascadia Code', monospace;
    font-size: 9pt;
    padding: 4px;
    background: transparent;
}
"""


