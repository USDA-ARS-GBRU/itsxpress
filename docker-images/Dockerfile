# Build directly on top of the official QIIME 2 Amplicon environment 
FROM quay.io/qiime2/qiime2-workshop:2026.7

LABEL org.opencontainers.image.source="https://github.com/usda-ars-gbru/itsxpress"

# Update/Install matching core dependencies within the container context
RUN conda install -y -c bioconda vsearch=2.30.4 hmmer=3.4 biopython=1.84 && \
    pip install pyzstd pytest

# Copy your local, pushed repository files cleanly into the working directory
WORKDIR /app
COPY . /app

# Install itsxpress directly into the active base environment
RUN pip install .

# Pre-refresh the plugin cache so the itsxpress plugin is instantly available to users
RUN qiime dev refresh-cache

# Set the default command to run itsxpress standalone CLI
CMD ["itsxpress"]
