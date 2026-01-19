# 🧬 FlowSOM Analyzer

**Application Desktop GUI pour l'analyse FlowSOM automatisée et interactive**

Application Python complète permettant de réaliser des analyses FlowSOM sur des données de cytométrie en flux (Flow Cytometry), avec comparaison entre échantillons sains et pathologiques.

![Python](https://img.shields.io/badge/Python-3.9+-blue)
![PyQt5](https://img.shields.io/badge/PyQt5-5.15+-green)
![FlowSOM](https://img.shields.io/badge/FlowSOM-saeyslab-orange)

---

## 📋 Fonctionnalités

### Interface Graphique (PyQt5)
- ✅ Design professionnel avec thème sombre moderne
- ✅ Chargement de dossiers FCS (Sain / Pathologique)
- ✅ Panneau de paramètres FlowSOM complet
- ✅ Visualisation intégrée Matplotlib avec barre d'outils

### Traitement FlowSOM
- ✅ Clustering SOM (Self-Organizing Map)
- ✅ Métaclustering par consensus
- ✅ Thread séparé (QThread) pour éviter le gel de l'interface
- ✅ Tracking de l'origine des fichiers (condition)

### Visualisations
- 📊 **Star Chart (MST)** : Vue Minimum Spanning Tree classique
- 📐 **Grid View** : Vue en grille SOM
- 🔥 **Heatmap** : Expression Z-score par métacluster
- 📈 **Distribution** : Comparaison Sain vs Pathologique
- 🎨 **Marker Expression** : Expression d'un marqueur spécifique

### Auto-Clustering
- 🔍 Détermination automatique du nombre optimal de métaclusters
- 📏 Basé sur le score de silhouette

### Export
- 💾 **FCS** : Export au format FCS standard avec colonnes cluster_id
- 📄 **CSV** : Export complet avec métadonnées
- 🖼️ **Image** : Export des figures (PNG, PDF, SVG)

---

## 🚀 Installation

### 1. Cloner ou télécharger le projet

```bash
cd FlowSom
```

### 2. Créer un environnement virtuel (recommandé)

```bash
python -m venv venv
# Windows
venv\Scripts\activate
# Linux/macOS
source venv/bin/activate
```

### 3. Installer les dépendances

```bash
pip install -r requirements.txt
```

### 4. Lancer l'application

```bash
python main.py
```

---

## 📖 Utilisation

### Étape 1 : Charger les données

1. Cliquez sur **"📁 Dossier Sain"** pour sélectionner le dossier contenant vos fichiers FCS sains
2. Cliquez sur **"📁 Dossier Pathologique"** pour sélectionner le dossier contenant vos fichiers FCS pathologiques

> **Note** : Vous pouvez charger un seul groupe si nécessaire

### Étape 2 : Configurer les paramètres

| Paramètre | Description | Valeur par défaut |
|-----------|-------------|-------------------|
| X Dimension | Largeur de la grille SOM | 10 |
| Y Dimension | Hauteur de la grille SOM | 10 |
| Nombre de Métaclusters | Nombre de groupes finaux | 10 |
| Seed | Graine pour reproductibilité | 42 |
| Auto-clustering | Détection automatique du nombre optimal | Désactivé |
| Exclure FSC/SSC/Time | Exclure les paramètres de scatter | Activé |

### Étape 3 : Lancer l'analyse

Cliquez sur **"▶️ Lancer l'Analyse FlowSOM"**

L'analyse s'exécute dans un thread séparé, vous pouvez suivre la progression en temps réel.

### Étape 4 : Explorer les résultats

Utilisez le menu déroulant pour basculer entre les différentes visualisations :
- **Star Chart (MST)** : Vue classique FlowSOM
- **Grid View** : Grille SOM avec métaclusters
- **Heatmap** : Profils d'expression
- **Distribution par Condition** : Comparaison statistique
- **Marker Expression** : Expression d'un marqueur spécifique

### Étape 5 : Exporter les résultats

- **📤 Exporter en FCS** : Fichier FCS avec colonnes de clustering ajoutées
- **📤 Exporter en CSV** : Tableau complet avec toutes les métadonnées
- **🖼️ Exporter la Figure** : Sauvegarder la visualisation courante

---

## 📦 Dépendances

```
PyQt5>=5.15.0
matplotlib>=3.5.0
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
scikit-learn>=1.0.0
anndata>=0.8.0
scanpy>=1.9.0
flowsom>=0.2.0
fcswrite>=0.6.0
flowio>=1.0.0
pytometry>=0.1.5 (optionnel)
```

---

## 🔧 Architecture du Code

```
FlowSom/
├── main.py              # Application complète
├── requirements.txt     # Dépendances Python
└── README.md            # Documentation
```

### Classes principales

| Classe | Description |
|--------|-------------|
| `FlowSOMApp` | Fenêtre principale PyQt5 |
| `FlowSOMWorker` | Thread de calcul FlowSOM |
| `MatplotlibCanvas` | Canvas intégré Matplotlib |

---

## 📚 Références

- **FlowSOM Python** : https://github.com/saeyslab/FlowSOM_Python
- **Documentation FlowSOM** : https://flowsom.readthedocs.io
- **Article original** : Van Gassen et al., "FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data", Cytometry Part A, 2015

---

## 📝 Licence

Ce projet est fourni à des fins éducatives et de recherche.

---

## 🤝 Contact

Pour toute question ou suggestion, veuillez ouvrir une issue sur le dépôt du projet.
