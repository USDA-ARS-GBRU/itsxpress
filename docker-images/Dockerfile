FROM condaforge/mambaforge:24.7.1-0

# # Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

LABEL org.opencontainers.image.source="https://github.com/usda-ars-gbru/itsxpress"
# Install conda dependencies
RUN conda update mamba -c conda-forge
RUN mamba install -c conda-forge python>=3.9
RUN mamba install -c bioconda vsearch=2.22.1 hmmer=3.1b2

# Copy the itsxpress package files and install dependencies
COPY .. /app
RUN cd /app && git clone https://github.com/USDA-ARS-GBRU/itsxpress.git && cd itsxpress && pip install .

# Set the default command to run itsxpress
CMD ["itsxpress"]
