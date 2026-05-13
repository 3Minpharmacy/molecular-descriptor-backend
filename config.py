"""
app/config.py
=============
Centralised configuration management using Pydantic BaseSettings.

All tuneable values are read from environment variables (with sensible
defaults for local development). Secrets should **never** be hard-coded;
use a .env file locally and platform-managed secrets in production.
"""

from __future__ import annotations

from typing import List

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Application-wide configuration resolved from the environment."""

    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        case_sensitive=False,
    )

    # ------------------------------------------------------------------
    # Application metadata
    # ------------------------------------------------------------------
    PROJECT_NAME: str = "MolecularMind — Computational Molecular Evaluation Engine"
    VERSION: str = "1.0.0"
    DESCRIPTION: str = (
        "Production-grade molecular descriptor engine for medicinal chemistry "
        "and AI-driven drug discovery."
    )

    # ------------------------------------------------------------------
    # Server
    # ------------------------------------------------------------------
    HOST: str = "0.0.0.0"
    PORT: int = 8000
    RELOAD: bool = False  # True only for local hot-reload development

    # ------------------------------------------------------------------
    # CORS
    # Comma-separated list in env: ALLOWED_ORIGINS="http://localhost:3000,https://myapp.com"
    # ------------------------------------------------------------------
    ALLOWED_ORIGINS: List[str] = ["*"]

    # ------------------------------------------------------------------
    # Logging
    # ------------------------------------------------------------------
    LOG_LEVEL: str = "INFO"

    # ------------------------------------------------------------------
    # Chemistry engine flags
    # ------------------------------------------------------------------
    # Maximum number of atoms allowed in a submitted molecule.
    # Prevents abuse with pathologically large SMILES strings.
    MAX_ATOM_COUNT: int = 500

    # Whether to attempt Mordred descriptor expansion (requires mordred install)
    ENABLE_MORDRED: bool = False


# Singleton — imported everywhere as `from app.config import settings`
settings = Settings()
