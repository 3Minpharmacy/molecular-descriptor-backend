"""
MolecularMind — Computational Molecular Evaluation Engine
==========================================================
Entry point for the FastAPI application.

This module initialises the application, registers middleware,
mounts routers, and configures lifecycle hooks. All domain
logic is deliberately kept out of this file; it is purely an
infrastructure concern.
"""

import logging
import time
from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from app.config import settings
from app.routes import descriptors

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.getLevelName(settings.LOG_LEVEL),
    format="%(asctime)s | %(levelname)-8s | %(name)s — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Application lifecycle
# ---------------------------------------------------------------------------

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup / shutdown hooks for resource management."""
    logger.info("MolecularMind engine starting — validating RDKit availability …")
    try:
        from rdkit import Chem  # noqa: F401
        logger.info("RDKit loaded successfully.")
    except ImportError as exc:
        logger.critical("RDKit is not installed: %s", exc)
        raise RuntimeError("RDKit is required but not available.") from exc
    yield
    logger.info("MolecularMind engine shutting down.")


# ---------------------------------------------------------------------------
# Application factory
# ---------------------------------------------------------------------------

def create_application() -> FastAPI:
    """Construct and configure the FastAPI application instance."""

    app = FastAPI(
        title=settings.PROJECT_NAME,
        description=(
            "A production-grade computational molecular evaluation engine for "
            "medicinal chemistry and AI-driven drug discovery workflows. "
            "Accepts SMILES strings and returns richly annotated descriptor "
            "profiles with medicinal chemistry reasoning."
        ),
        version=settings.VERSION,
        docs_url="/docs",
        redoc_url="/redoc",
        openapi_url="/openapi.json",
        lifespan=lifespan,
    )

    # -----------------------------------------------------------------------
    # CORS — allow configured origins (or wildcard in development)
    # -----------------------------------------------------------------------
    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.ALLOWED_ORIGINS,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # -----------------------------------------------------------------------
    # Request latency logging middleware
    # -----------------------------------------------------------------------
    @app.middleware("http")
    async def log_requests(request: Request, call_next):
        start = time.perf_counter()
        response = await call_next(request)
        elapsed_ms = (time.perf_counter() - start) * 1_000
        logger.info(
            "%s %s → %d (%.1f ms)",
            request.method,
            request.url.path,
            response.status_code,
            elapsed_ms,
        )
        return response

    # -----------------------------------------------------------------------
    # Global exception handler — catch unhandled errors gracefully
    # -----------------------------------------------------------------------
    @app.exception_handler(Exception)
    async def unhandled_exception_handler(request: Request, exc: Exception):
        logger.exception("Unhandled exception on %s: %s", request.url.path, exc)
        return JSONResponse(
            status_code=500,
            content={
                "detail": "An unexpected internal error occurred.",
                "path": str(request.url.path),
            },
        )

    # -----------------------------------------------------------------------
    # Routers
    # -----------------------------------------------------------------------
    app.include_router(
        descriptors.router,
        prefix="/api/v1",
        tags=["Molecular Descriptors"],
    )

    # -----------------------------------------------------------------------
    # Health probe — used by Render / Railway / Kubernetes liveness checks
    # -----------------------------------------------------------------------
    @app.get("/health", tags=["Health"], summary="Liveness probe")
    async def health_check():
        return {"status": "ok", "service": settings.PROJECT_NAME, "version": settings.VERSION}

    return app


app = create_application()
