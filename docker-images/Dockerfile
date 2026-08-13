# Build directly on top of the official QIIME 2 environment
FROM quay.io/qiime2/qiime2:2026.7

LABEL org.opencontainers.image.source="https://github.com"

# Update/Install ONLY the missing additions within the active container context
# vsearch, hmmer, and pytest are skipped because the base image already includes them
RUN conda install -y -c conda-forge pyzstd==0.19.1 && \
    pip install biopython==1.88

# Copy your local, pushed repository files cleanly into the working directory
WORKDIR /app
COPY . /app

# Ensure build dependencies are ready and install itsxpress without build isolation
# This prevents setuptools_scm from failing due to missing build isolation requirements
RUN pip install --upgrade pip setuptools wheel setuptools_scm && \
    pip install --no-build-isolation .

# Pre-refresh the plugin cache so the itsxpress plugin is instantly available to users
RUN qiime dev refresh-cache

# Set the default command to run itsxpress standalone CLI
CMD ["itsxpress"]
