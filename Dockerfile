FROM continuumio/miniconda3:4.10.3

# Install system dependencies
RUN apt-get update && \
    apt-get install -y build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

LABEL org.opencontainers.image.source="https://github.com/usda-ars-gbru/itsxpress"
# Install conda dependencies
RUN conda install -c bioconda vsearch=2.22.1 hmmer=3.1b2

# Copy the itsxpress package files and install dependencies
COPY . /app
WORKDIR /app
RUN pip install --no-cache-dir .

# Set the default command to run itsxpress
CMD ["itsxpress"]