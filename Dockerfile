# ===========================================================
# echolocatoR Docker Image
# ===========================================================
# Automated genomic fine-mapping pipeline with all dependencies
# pre-installed: R packages, Python/conda env, and system tools.
#
# Build:
#   docker build -t echolocator .
# Run:
#   docker run -d -e ROOT=true -e PASSWORD=bioc \
#     -v ~/Desktop:/Desktop -p 8788:8787 echolocator
# Then open http://localhost:8788 (user: rstudio, password: bioc)
# ===========================================================

ARG BASE_IMAGE=ghcr.io/bioconductor/bioconductor_docker:devel
FROM ${BASE_IMAGE}

# ---- System dependencies ----
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git-core \
    libcurl4-openssl-dev \
    libgit2-dev \
    libicu-dev \
    libssl-dev \
    make \
    pandoc \
    zlib1g-dev \
    xfonts-100dpi \
    xfonts-75dpi \
    biber \
    libsbml5-dev \
    qpdf \
    cmake \
    g++ \
    wget \
    curl \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ---- Install Miniforge (conda) ----
ENV CONDA_DIR=/opt/miniforge
RUN curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    -o /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniforge.sh
ENV PATH="${CONDA_DIR}/bin:${PATH}"

# ---- Create echoR_mini conda environment ----
COPY inst/conda/echoR_mini.yml /tmp/echoR_mini.yml
# If the yml doesn't exist in this repo, fetch from echoconda
RUN if [ ! -f /tmp/echoR_mini.yml ]; then \
      curl -fsSL https://raw.githubusercontent.com/RajLabMSSM/echoconda/main/inst/conda/echoR_mini.yml \
        -o /tmp/echoR_mini.yml; \
    fi && \
    conda env create -f /tmp/echoR_mini.yml && \
    conda clean -afy && \
    rm /tmp/echoR_mini.yml

# Make conda available to R
ENV RETICULATE_MINICONDA_PATH="${CONDA_DIR}"
RUN echo "RETICULATE_MINICONDA_PATH=${CONDA_DIR}" >> /usr/local/lib/R/etc/Renviron.site && \
    echo "PATH=${CONDA_DIR}/bin:${CONDA_DIR}/envs/echoR_mini/bin:\${PATH}" >> /usr/local/lib/R/etc/Renviron.site

# ---- Install echolocatoR + all R dependencies (Imports + Suggests) ----
ARG PKG=echolocatoR
RUN echo ${PKG}
RUN mkdir -p /${PKG}
ADD . /${PKG}
WORKDIR /${PKG}

RUN --mount=type=secret,id=GITHUB_TOKEN,env=GITHUB_PAT \
    Rscript -e ' \
    options(download.file.method="libcurl", crayon.enabled=TRUE, timeout=2000); \
    if(!require("BiocManager")) install.packages("BiocManager"); \
    if(!require("remotes")) install.packages("remotes"); \
    repos <- BiocManager::repositories(); \
    remotes::install_local(repos=repos, dependencies=TRUE, \
                           build_vignettes=FALSE, upgrade=TRUE, force=TRUE);'

# ---- Install genetics.binaRies (GCTA + plink + bcftools) ----
RUN Rscript -e ' \
    if(!require("genetics.binaRies", quietly=TRUE)) \
      remotes::install_github("RajLabMSSM/genetics.binaRies"); \
    message("GCTA: ", genetics.binaRies::get_gcta_binary()); \
    message("plink: ", genetics.binaRies::get_plink_binary()); \
    message("bcftools: ", genetics.binaRies::get_bcftools_binary());'

# ---- Download FINEMAP binary ----
RUN Rscript -e ' \
    fp <- echofinemap:::FINEMAP_find_executable(verbose=TRUE); \
    message("FINEMAP installed at: ", fp);'

# ---- Compile PAINTOR from submodule ----
RUN Rscript -e ' \
    tryCatch({ \
      echofinemap:::PAINTOR_install(verbose=TRUE); \
      message("PAINTOR installed successfully."); \
    }, error = function(e) message("PAINTOR install failed: ", e$message));'

# ---- Verify installation ----
RUN Rscript -e ' \
    library(echolocatoR); \
    message("echolocatoR ", packageVersion("echolocatoR"), " loaded successfully"); \
    message("Conda envs: ", paste(echoconda::list_envs(), collapse=", ")); \
    '

# ---- Clean up build directory ----
RUN rm -rf /${PKG}
WORKDIR /home/rstudio
