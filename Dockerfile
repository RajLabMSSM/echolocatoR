# ===========================================================
# echolocatoR Docker Image
# ===========================================================
# Automated genomic fine-mapping pipeline with all dependencies
# pre-installed: R packages, Python/conda env, and system tools.
#
# Uses echofinemap container as base — it already includes
# Bioconductor, echofinemap, echodata, echotabix, echoLD,
# echoconda, and all their dependencies.
#
# Build (local):
#   docker build --build-arg GITHUB_PAT=<token> -t echolocator .
# Run:
#   docker run -d -e ROOT=true -e PASSWORD=bioc \
#     -v ~/Desktop:/Desktop -p 8788:8787 echolocator
# Then open http://localhost:8788 (user: rstudio, password: bioc)
# ===========================================================

ARG BASE_IMAGE=ghcr.io/rajlabmssm/echofinemap:latest
FROM ${BASE_IMAGE}

# ---- System dependencies (extras not in base) ----
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    cmake \
    g++ \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ---- Install Miniforge (conda) ----
ENV CONDA_DIR=/opt/miniforge
RUN ARCH=$(uname -m) && \
    curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${ARCH}.sh" \
    -o /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniforge.sh
ENV PATH="${CONDA_DIR}/bin:${PATH}"

# ---- Create echoR_mini conda environment ----
COPY inst/conda/echoR_mini.yml /tmp/echoR_mini.yml
RUN if [ ! -f /tmp/echoR_mini.yml ] || [ ! -s /tmp/echoR_mini.yml ]; then \
      curl -fsSL https://raw.githubusercontent.com/RajLabMSSM/echoconda/main/inst/conda/echoR_mini.yml \
        -o /tmp/echoR_mini.yml; \
    fi && \
    conda env create -f /tmp/echoR_mini.yml && \
    conda clean -afy && \
    rm /tmp/echoR_mini.yml

# Make conda available to R and system-wide
ENV RETICULATE_MINICONDA_PATH="${CONDA_DIR}"
ENV CONDA_DEFAULT_ENV="echoR_mini"
RUN echo "RETICULATE_MINICONDA_PATH=${CONDA_DIR}" >> /usr/local/lib/R/etc/Renviron.site && \
    echo "PATH=${CONDA_DIR}/bin:${CONDA_DIR}/envs/echoR_mini/bin:\${PATH}" >> /usr/local/lib/R/etc/Renviron.site && \
    ln -sf ${CONDA_DIR}/bin/conda /usr/local/bin/conda

# ---- Install echolocatoR + remaining R deps ----
ARG PKG=echolocatoR
RUN mkdir -p /${PKG}
ADD . /${PKG}
WORKDIR /${PKG}

# On GHA, rworkflows passes GITHUB_TOKEN via --mount=type=secret.
# For local builds, pass --build-arg GITHUB_PAT=<token>.
ARG GITHUB_PAT=""
ENV GITHUB_PAT=${GITHUB_PAT}
RUN --mount=type=secret,id=GITHUB_TOKEN \
    export GITHUB_PAT="${GITHUB_PAT:-$(cat /run/secrets/GITHUB_TOKEN 2>/dev/null || echo '')}" && \
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
      remotes::install_github("MRCIEU/genetics.binaRies"); \
    message("GCTA: ", genetics.binaRies::get_gcta_binary()); \
    message("plink: ", genetics.binaRies::get_plink_binary()); \
    message("bcftools: ", genetics.binaRies::get_bcftools_binary());'

# ---- Download FINEMAP binary ----
RUN Rscript -e ' \
    fp <- echofinemap:::FINEMAP_find_executable(verbose=TRUE); \
    message("FINEMAP installed at: ", fp);'

# ---- Compile PAINTOR from echofinemap submodule ----
RUN Rscript -e ' \
    paintor_src <- system.file("tools", "PAINTOR_V3.0", package="echofinemap"); \
    if(nchar(paintor_src) > 0 && file.exists(file.path(paintor_src, "install.sh"))) { \
      message("Compiling PAINTOR from: ", paintor_src); \
      setwd(paintor_src); \
      system("bash install.sh"); \
      message("PAINTOR installed successfully."); \
    } else { \
      message("NOTE: PAINTOR source not bundled in installed echofinemap."); \
      message("Users can compile later: echofinemap:::PAINTOR_install()"); \
    }'

# ---- Verify installation ----
RUN Rscript -e ' \
    library(echolocatoR); \
    message("echolocatoR ", packageVersion("echolocatoR"), " loaded successfully"); \
    message("FINEMAP: ", tryCatch(echofinemap:::FINEMAP_find_executable(verbose=FALSE), error=function(e) "not found")); \
    message("Conda: ", tryCatch(echoconda::find_conda(verbose=FALSE), error=function(e) "not found")); \
    message("echoR_mini env: ", tryCatch(echoconda::env_exists("echoR_mini"), error=function(e) "not found")); \
    '

# ---- Clean up ----
RUN rm -rf /${PKG}
WORKDIR /home/rstudio
