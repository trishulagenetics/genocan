FROM nfcore/base
MAINTAINER Andries van Tonder <trishulagenetics@gmail.com>
LABEL authors="trishulagenetics@gmail.com" \
    description="Docker image containing all requirements for trishulagenetics/genocan pipeline"

COPY environment.yml /
RUN conda env update -n root -f /environment.yml && conda clean -a