From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Andries van Tonder <trishulagenetics@gmail.com>
    DESCRIPTION Singularity image containing all requirements for trishulagenetics/genocan pipeline
    VERSION 0.1dev

%files
    environment.yml /

%post
    /opt/conda/bin/conda env update -n root -f /environment.yml
    /opt/conda/bin/conda clean -a