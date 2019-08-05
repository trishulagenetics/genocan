FROM debian:jessie

MAINTAINER Andries van Tonder <trishulagenetics@gmail.com>
LABEL authors="trishulagenetics@gmail.com" \
    description="Docker image containing all requirements for trishulagenetics/genocan pipeline"

RUN apt-get update && apt-get install --yes --no-install-recommends \
	build-essential \
	wget \
	unzip \
	perl \
	cmake \
	graphviz \
	python-dev \
	python-pip \
	python3-dev \
	python3-pip \
	cpanminus \
	libfreetype6-dev \
	pkg-config \
	libpng12-dev \
	libfindbin-libs-perl \
	python3-matplotlib \
	default-jre \
	r-base \
	libncurses-dev \
	libbz2-dev \
	liblzma-dev \
	zlib1g-dev \
	tabix \
	&& rm -rf /var/lib/apt/lists/*

RUN cpanm --force CPAN::Meta \
	#Path::Tiny \
	#File::Copy::Recursive::Reduced \
	FindBin::Real

RUN echo "deb http://cran.rstudio.com/bin/linux/debian jessie-cran3/" >>  /etc/apt/sources.list \
 && apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480 \
 && apt-get update --fix-missing \
 && apt-get -y install r-base

WORKDIR /tools
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
RUN unzip Trimmomatic-0.36
RUN mv Trimmomatic-0.36 trimmomatic

ENTRYPOINT ["java", "-jar","/tools/trimmomatic/trimmomatic-0.36.jar"]

ENV DST=/tmp

# FASTQC
ENV URL=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
ENV ZIP=fastqc_v0.11.8.zip

RUN wget $URL/$ZIP -O $DST/$ZIP && \
  unzip - $DST/$ZIP -d $DST && \
  rm $DST/$ZIP && \
  cd $DST/FastQC && \
  chmod 755 fastqc && \
  ln -s $DST/FastQC/fastqc /usr/local/bin/fastqc
  
ENV PATH /usr/local/bin:$PATH

CMD fastqc

## Install MultiQC:
RUN pip3 install 'networkx==2.2'

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
RUN pip3 install multiqc

CMD multiqc

RUN R -e 'install.packages("markdown", dependencies=TRUE, repos="http://cloud.r-project.org/")'

ENV BWA_SOURCE="bwa-0.7.15.tar.bz2" \
    BWA_VERSION="0.7.15" \
    BWA_BIN="bwa" \
    BWA_DEST="/usr/local/bin/bwa"

RUN wget https://github.com/lh3/bwa/releases/download/v$BWA_VERSION/$BWA_SOURCE -O /opt/$BWA_SOURCE \
    && tar -xvf /opt/$BWA_SOURCE -C /opt \
    && cd /opt/bwa-$BWA_VERSION \
    && make \
    && ln -s /opt/bwa-$BWA_VERSION/$BWA_BIN $BWA_DEST \
    && rm /opt/$BWA_SOURCE

ENV SAMTOOLS_SOURCE="samtools-1.4.1.tar.bz2" \
    SAMTOOLS_VERSION="1.4.1" \
    SAMTOOLS_BIN="samtools" \
    SAMTOOLS_DEST="/usr/local/bin/samtools"

RUN wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/$SAMTOOLS_SOURCE -O /opt/$SAMTOOLS_SOURCE \
    && tar -xvf /opt/$SAMTOOLS_SOURCE -C /opt \
    && cd /opt/samtools-$SAMTOOLS_VERSION \
    && ./configure \
    && make \
    && ln -s /opt/samtools-$SAMTOOLS_VERSION/$SAMTOOLS_BIN $SAMTOOLS_DEST \
    && rm /opt/$SAMTOOLS_SOURCE

ENV BCFTOOLS_SOURCE="bcftools-1.4.1.tar.bz2" \
    BCFTOOLS_VERSION="1.4.1" \
    BCFTOOLS_BIN="bcftools" \
    BCFTOOLS_DEST="/usr/local/bin/bcftools"

RUN wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/$BCFTOOLS_SOURCE -O /opt/$BCFTOOLS_SOURCE \
    && tar -xvf /opt/$BCFTOOLS_SOURCE -C /opt \
    && cd /opt/bcftools-$BCFTOOLS_VERSION \
    && make \
    && ln -s /opt/bcftools-$BCFTOOLS_VERSION/$BCFTOOLS_BIN $BCFTOOLS_DEST \
    && rm /opt/$BCFTOOLS_SOURCE

ENV ZIP=vcftools-0.1.15.tar.gz
ENV URL=https://github.com/vcftools/vcftools/releases/download/v0.1.15/
ENV FOLDER=vcftools-0.1.15
ENV DST=/tmp

RUN wget $URL/$ZIP -O $DST/$ZIP && \
  tar xvf $DST/$ZIP -C $DST && \
  rm $DST/$ZIP && \
  cd $DST/$FOLDER && \
  ./configure && \
  make && \
  make install && \
  cd / && \
  rm -rf $DST/$FOLDER