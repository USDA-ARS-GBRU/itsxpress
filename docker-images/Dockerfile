# Build directly on top of the official QIIME 2 environment
FROM quay.io/qiime2/qiime2:2026.7

LABEL org.opencontainers.image.source="https://github.com"

# 1. Update/Install ONLY the missing additions directly within the active QIIME environment
RUN conda run -n qiime2-core conda install -y -c conda-forge pyzstd==0.19.1 && \
    conda run -n qiime2-core pip install biopython==1.88

# Copy your local repository files cleanly into the working directory
WORKDIR /app
COPY . /app

# 2. Upgrade build tools inside the QIIME environment
RUN conda run -n qiime2-core pip install --upgrade pip setuptools wheel setuptools_scm

# 3. Install itsxpress without build isolation to prevent setuptools_scm lookup failures
RUN conda run -n qiime2-core pip install --no-build-isolation .

# 4. Pre-refresh the plugin cache inside the target environment
RUN conda run -n qiime2-core qiime dev refresh-cache

# Set the default command to run itsxpress standalone CLI
CMD ["itsxpress"]
