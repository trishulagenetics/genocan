FROM debian:jessie

MAINTAINER Andries van Tonder <trishulagenetics@gmail.com>
LABEL authors="trishulagenetics@gmail.com" \
    description="Docker image containing all requirements for trishulagenetics/genocan pipeline"

RUN apt-get update && apt-get install --yes --no-install-recommends \
	wget \
	wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
	RUN unzip Trimmomatic-0.36
	RUN mv Trimmomatic-0.36 trimmomatic


ENV DST=/tmp
ENV URL=https://github.com/agordon/libgtextutils/releases/download/0.7/

RUN wget $URL/libgtextutils-0.7.tar.gz -O $DST/libgtextutils-0.7.tar.gz && \
  tar -xvf $DST/libgtextutils-0.7.tar.gz -C $DST && \
  rm $DST/libgtextutils-0.7.tar.gz && \
  cd $DST/libgtextutils-0.7 && \
  ./configure && \
  make && \
  make install && \
  cd / && \
rm -rf $DST/libgtextutils-0.7

ENV URL=https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/

# FastXtools
RUN wget $URL/fastx_toolkit-0.0.14.tar.bz2 -O $DST/fastx_toolkit-0.0.14.tar.bz2 && \
  tar -xvf $DST/fastx_toolkit-0.0.14.tar.bz2 -C $DST && \
  rm $DST/fastx_toolkit-0.0.14.tar.bz2 && \
  cd $DST/fastx_toolkit-0.0.14 && \
  ./configure && \
  make && \
  make install && \
  cd / && \
  rm -rf $DST/fastx_toolkit-0.0.14

RUN ldconfig 

# FASTQC
ENV URL=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
ENV ZIP=fastqc_v0.11.5.zip

RUN wget $URL/$ZIP -O $DST/$ZIP && \
  unzip - $DST/$ZIP -d $DST && \
  rm $DST/$ZIP && \
  cd $DST/FastQC && \
  chmod 755 fastqc && \
  ln -s $DST/FastQC/fastqc /usr/local/bin/fastqc
  
ENV PATH /usr/local/bin:$PATH

CMD fastqc

RUN apk add --no-cache bash

# Add the MultiQC source files to the container
ADD . /usr/src/multiqc
WORKDIR /usr/src/multiqc

# Remove matplotlib requirement needed for py2 support
# TODO: We can get rid of this when MultiQC is py3 only
RUN sed -i 's/matplotlib>=2.1.1,<3.0.0/matplotlib>=2.1.1/g' setup.py

# Install MultiQC
RUN python setup.py install

# Set up entrypoint and cmd for easy docker usage
ENTRYPOINT [ "multiqc" ]
CMD multiqc

WORKDIR /data/
