# Dockerfile for package buliding and running this package
FROM rocker/r-ver:3.6.1

RUN apt-get update && apt-get -y install \
  texlive-latex-base \
  texlive-fonts-extra \
  texinfo \
  libxml2-dev

RUN R -e "install.packages('remotes')"

COPY DESCRIPTION /tmp/DESCRIPTION

RUN R -e "library('remotes'); install_deps('/tmp', dependencies = TRUE)"
