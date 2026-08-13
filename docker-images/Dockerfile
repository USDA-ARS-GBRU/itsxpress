# Build directly on top of the official QIIME 2 environment
FROM quay.io/qiime2/qiime2:2026.7

LABEL org.opencontainers.image.source="https://github.com"

# Update/Install ONLY the missing additions within the active container context
# vsearch, hmmer, and pytest are skipped because the base image already includes them
RUN conda install -y -c conda-forge pyzstd==0.19.1 && \
    python -m pip install biopython==1.88

# Copy your local, pushed repository files cleanly into the working directory
WORKDIR /app
COPY . /app

# Accept the version dynamically from GitHub Actions (Defaults to 2.0.0 if run manually)
ARG DOCKER_IMAGE_VERSION="2.0.0"

# Best Practice: Standard pip install. The inline environment variable instructs 
# setuptools_scm to skip git discovery and build the full wheel natively.
RUN SETUPTOOLS_SCM_PRETEND_VERSION=$DOCKER_IMAGE_VERSION python -m pip install .

# Pre-refresh the plugin cache so the itsxpress plugin is instantly available to users
RUN qiime dev refresh-cache

# Set the default command to run itsxpress standalone CLI
CMD ["itsxpress"]
