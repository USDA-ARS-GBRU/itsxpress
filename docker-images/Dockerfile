# Build directly on top of the official QIIME 2 environment
FROM quay.io/qiime2/qiime2:2026.7

LABEL org.opencontainers.image.source="https://github.com"

# Update/Install ONLY the missing additions within the active container context
RUN conda install -y -c conda-forge pyzstd==0.19.1 && \
    python -m pip install biopython==1.88

# Copy your local, pushed repository files cleanly into the working directory
WORKDIR /app
COPY . /app

# Accept the version dynamically from GitHub Actions
ARG DOCKER_IMAGE_VERSION="0.0.0"

# Best Practice Fix: Sanitize the version variable. If it starts with 'sha-', 
# format it into a PEP 440 compliant local version identifier '0.0.0+sha.edece1c2'.
RUN CLEAN_VERSION=$(echo "$DOCKER_IMAGE_VERSION" | sed 's/^sha-/0.0.0+sha./') && \
    SETUPTOOLS_SCM_PRETEND_VERSION="$CLEAN_VERSION" python -m pip install .

# Pre-refresh the plugin cache so the itsxpress plugin is instantly available to users
RUN qiime dev refresh-cache

# Set the default command to run itsxpress standalone CLI
CMD ["itsxpress"]
