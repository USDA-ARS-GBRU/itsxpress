FROM continuumio/miniconda3:4.10.3

# Install system dependencies
RUN apt-get update && \
    apt-get install -y build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install conda dependencies
RUN conda install -c bioconda vsearch hmmer

# Copy the itsxpress package files and install dependencies
COPY . /app
WORKDIR /app
RUN pip install --no-cache-dir .

# Set the default command to run itsxpress
CMD ["itsxpress"]