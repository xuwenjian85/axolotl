#
# NOTE: THIS DOCKERFILE IS GENERATED VIA "update.sh"
#
# PLEASE DO NOT EDIT IT DIRECTLY.
#

ARG BASE_IMAGE=rstudio/r-base
FROM ${BASE_IMAGE}:jammy

ARG R_VERSION=4.2.3
ARG OS_IDENTIFIER=ubuntu-2204

# Install R
RUN wget https://cdn.posit.co/r/${OS_IDENTIFIER}/pkgs/r-${R_VERSION}_1_amd64.deb && \
    apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -f -y ./r-${R_VERSION}_1_amd64.deb && \
    ln -s /opt/R/${R_VERSION}/bin/R /usr/bin/R && \
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/bin/Rscript && \
    ln -s /opt/R/${R_VERSION}/lib/R /usr/lib/R && \
    rm r-${R_VERSION}_1_amd64.deb 
    # && \
    # rm -rf /var/lib/apt/lists/*

RUN apt install -y libssl-dev
RUN apt-get install -y libpng-dev libxml2-dev

RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
RUN R -e "install.packages('httr',dependencies=TRUE, repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
RUN R -e "install.packages('png',dependencies=TRUE, repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"

RUN R -e "\
options('repos' = c(CRAN='https://mirrors.tuna.tsinghua.edu.cn/CRAN/'));\
options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');\
BiocManager::install('OUTRIDER',update = FALSE)"

# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# BiocManager::install(c("KEGGREST","plotly"))
# BiocManager::install(c("AnnotationDbi"))
# BiocManager::install("OUTRIDER",update = FALSE)


# AnnotationDbi
# CMD ["R"]

