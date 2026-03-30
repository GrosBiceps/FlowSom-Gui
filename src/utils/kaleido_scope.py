"""kaleido_scope.py — Scope kaleido persistant pour les exports Plotly → image.

Kaleido lance par défaut un processus Chromium à chaque appel de write_image()
ou pio.to_image(), ce qui coûte ~2-4 s de startup par figure.

Ce module initialise le scope kaleido UNE SEULE FOIS par processus et le réutilise
pour tous les appels suivants, réduisant le surcoût à ~0 après la première figure.

Utilisation :
    from flowsom_pipeline_pro.src.utils.kaleido_scope import ensure_kaleido_scope
    ensure_kaleido_scope()   # appel idempotent, à faire avant write_image()

    fig.write_image(path, format="jpg", ...)   # réutilise le scope existant
"""
from __future__ import annotations

import logging

_logger = logging.getLogger("utils.kaleido_scope")
_scope_initialized: bool = False


def ensure_kaleido_scope() -> bool:
    """Initialise le scope kaleido persistant si ce n'est pas déjà fait.

    Returns:
        True si kaleido est disponible et le scope est actif, False sinon.
    """
    global _scope_initialized
    if _scope_initialized:
        return True

    try:
        import plotly.io as pio

        scope = getattr(pio, "kaleido", None)
        if scope is None:
            _logger.debug("pio.kaleido non disponible — version kaleido trop ancienne.")
            return False

        # Forcer l'initialisation du scope (démarre le processus Chromium une fois)
        _ = scope.scope  # accès à l'attribut suffit à instancier le scope
        _scope_initialized = True
        _logger.debug("Scope kaleido persistant initialisé.")
        return True

    except Exception as exc:
        _logger.debug("Impossible d'initialiser le scope kaleido persistant: %s", exc)
        return False


def warm_up_kaleido() -> bool:
    """Rend une image vide pour forcer le démarrage de Chromium en avance.

    À appeler une seule fois au début de la pipeline, avant les vrais exports.
    Évite que le délai de startup soit imputé à la première vraie figure.

    Returns:
        True si le warm-up a réussi, False sinon.
    """
    try:
        import plotly.graph_objects as go
        import plotly.io as pio

        fig = go.Figure()
        pio.to_image(fig, format="png", width=10, height=10)
        _logger.debug("Kaleido warm-up effectué.")
        return True
    except Exception as exc:
        _logger.debug("Kaleido warm-up échoué (non bloquant): %s", exc)
        return False
