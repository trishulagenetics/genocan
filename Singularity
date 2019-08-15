From:trishulagenetics/genocan
Bootstrap:docker

%labels
    MAINTAINER Andries van Tonder <trishulagenetics@gmail.com>
    DESCRIPTION Singularity image containing all requirements for trishulagenetics/genocan pipeline
    VERSION 0.1dev

%environment
    PATH=/opt/conda/envs/trishulagenetics-genocan-0.1dev/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a