FROM nfcore/base

LABEL description="Docker image containing all requirements for trishulagenetics/genocan"
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/trishula-genocan-0.1dev/bin:$PATH