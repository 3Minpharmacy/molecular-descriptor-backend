# =============================================================================
# MolecularMind — Production Dockerfile
# =============================================================================
# Multi-stage build:
#   Stage 1 (builder) — install Python deps into an isolated venv
#   Stage 2 (runtime) — minimal image with only what's needed to run
#
# Build:  docker build -t molecularmind:latest .
# Run:    docker run -p 8000:8000 molecularmind:latest
# =============================================================================

# ---- Stage 1: dependency builder --------------------------------------------
FROM python:3.11-slim AS builder

# Install system-level build dependencies for RDKit and other compiled packages
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libboost-all-dev \
        libxrender1 \
        libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Create an isolated virtual environment to carry into the runtime stage
ENV VIRTUAL_ENV=/opt/venv
RUN python -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Upgrade pip inside the venv
RUN pip install --upgrade pip wheel

# Copy and install Python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt


# ---- Stage 2: minimal runtime image ----------------------------------------
FROM python:3.11-slim AS runtime

LABEL org.opencontainers.image.title="MolecularMind"
LABEL org.opencontainers.image.description="Computational Molecular Evaluation Engine"
LABEL org.opencontainers.image.version="1.0.0"

# Runtime system libraries required by RDKit
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxrender1 \
        libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Copy the pre-built virtual environment from the builder stage
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Create a non-root user for security best practice
RUN groupadd --gid 1001 appgroup && \
    useradd --uid 1001 --gid appgroup --shell /bin/bash --create-home appuser

WORKDIR /app

# Copy application source — as non-root
COPY --chown=appuser:appgroup app/ ./app/

USER appuser

# Expose the application port
EXPOSE 8000

# Health check — Docker / Kubernetes liveness probe
HEALTHCHECK --interval=30s --timeout=10s --start-period=15s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/health')"

# Production entry point — Uvicorn with Gunicorn-compatible process management
CMD ["uvicorn", "app.main:app", \
     "--host", "0.0.0.0", \
     "--port", "8000", \
     "--workers", "2", \
     "--log-level", "info", \
     "--access-log"]
