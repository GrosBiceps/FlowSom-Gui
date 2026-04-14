# -*- coding: utf-8 -*-
"""
exceptions.py — Exceptions cliniques et architecturales de FlowSOM Pipeline Pro.

Ces exceptions sont conçues pour être interceptées par l'UI (PyQt5) afin
d'afficher un message d'erreur explicite à l'utilisateur clinicien, sans
crash silencieux ni fallback dangereux.
"""

from __future__ import annotations


class ClinicalMathError(Exception):
    """
    Erreur mathématique à impact clinique potentiel.

    Levée lorsqu'une opération mathématique critique (inversion de matrice,
    calcul de covariance, distance de Mahalanobis) ne peut pas être effectuée
    de manière fiable. Le résultat de scoring serait invalide si l'on
    continuait silencieusement.

    Attributes:
        message: Description lisible par l'utilisateur clinicien.
        condition_number: Conditionnement de la matrice (si applicable).
            Un conditionnement > 1e12 indique une matrice quasi-singulière.
        details: Informations techniques supplémentaires pour le log.
    """

    def __init__(
        self,
        message: str,
        condition_number: float | None = None,
        details: str | None = None,
    ) -> None:
        self.message = message
        self.condition_number = condition_number
        self.details = details
        full_msg = message
        if condition_number is not None:
            full_msg += f" (cond={condition_number:.2e})"
        if details:
            full_msg += f" — {details}"
        super().__init__(full_msg)


class PanelConfigError(Exception):
    """
    Erreur de configuration du panel de marqueurs.

    Levée lorsqu'aucun panel n'est sélectionné ou que le fichier de
    configuration du panel est absent ou invalide. Empêche l'exécution
    avec des poids de scoring non définis.
    """

    def __init__(self, message: str, panel_path: str | None = None) -> None:
        self.message = message
        self.panel_path = panel_path
        full_msg = message
        if panel_path:
            full_msg += f" (chemin : {panel_path})"
        super().__init__(full_msg)


class PipelineStageError(Exception):
    """
    Erreur dans une étape du pipeline.

    Levée pour signaler qu'une étape spécifique du pipeline (gating,
    clustering, MRD) a échoué de manière non récupérable.
    """

    def __init__(self, stage: str, message: str) -> None:
        self.stage = stage
        self.message = message
        super().__init__(f"[{stage}] {message}")
