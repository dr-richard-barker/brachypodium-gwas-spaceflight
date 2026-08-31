FROM continuumio/miniconda3:latest

LABEL maintainer="Richard Barker <richard.barker@phylo.com>"
LABEL description="Brachypodium GWAS-Spaceflight Integration analysis environment"

WORKDIR /workspace

# Copy environment file and install dependencies
COPY environment.yml .
RUN conda env create -f environment.yml -n brachypodium-gwas-spaceflight && conda clean -afy

# Copy project files
COPY . /workspace/brachypodium-gwas-spaceflight

# Activate environment
ENV PATH="/opt/conda/envs/brachypodium-gwas-spaceflight/bin:$PATH"

# Set working directory to project root
WORKDIR /workspace/brachypodium-gwas-spaceflight

# Default command
CMD ["bash"]
