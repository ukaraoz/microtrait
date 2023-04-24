FROM rocker/r-ubuntu
MAINTAINER Ulas Karaoz <ukaraoz@lbl.gov>
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

ENV HMMER_VERSION='3.3.2'
ENV prodigal_VERSION='2.6.3'

RUN apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive apt-get -y install tzdata && \
    apt-get -y install make git ant curl gnupg-agent && \
    apt-get -y install apt-transport-https ca-certificates software-properties-common
#update-java-alternatives -s java-1.8.0-openjdk-amd64

RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add - && \
    add-apt-repository \
    "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
    $(lsb_release -cs) \
    stable" && \
    apt-get -y update

WORKDIR /tmp
RUN \
    curl http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz > hmmer-${HMMER_VERSION}.tar.gz && \
    tar xfz hmmer-${HMMER_VERSION}.tar.gz && \
    ln -s hmmer-${HMMER_VERSION} hmmer && \
    rm -f hmmer-${HMMER_VERSION}.tar.gz
WORKDIR /tmp/hmmer
RUN \
    ./configure --prefix /usr/local && \
    make && \
    make install

WORKDIR /tmp
RUN curl --location https://github.com/hyattpd/Prodigal/releases/download/v${prodigal_VERSION}/prodigal.linux > prodigal && \
    chmod +x prodigal && \
    cp prodigal /usr/local/bin/

RUN apt-get install -y --no-install-recommends \
    r-cran-dplyr r-cran-ape r-cran-assertthat r-cran-checkmate \
    r-cran-futile.logger r-cran-gtools r-cran-lazyeval r-cran-magrittr \
    r-cran-pheatmap r-cran-R.utils r-cran-stringr r-cran-tibble \
    r-cran-tidyr r-cran-readr r-cran-rcolorbrewer r-cran-corrplot \
    r-cran-doparallel r-cran-tictoc r-cran-biocmanager r-cran-devtools \
    r-bioc-biostrings r-bioc-complexheatmap r-cran-ggplot2 \
    r-cran-optparse
RUN R -e "install.packages(c('kmed'))"
#
RUN curl --location https://github.com/ukaraoz/microtrait/releases/download/kb/microtrait_1.0.0-kb.tar.gz > /tmp/microtrait_1.0.0-kb.tar.gz && \
    R CMD INSTALL /tmp/microtrait_1.0.0-kb.tar.gz

RUN mkdir -p /data
COPY ./testdata/* /data
COPY ./run_microtrait.R /usr/local/bin


