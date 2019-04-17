# This base image is meant for Java applications.
# It is currently set to use Java 8.

# https://registry.hub.docker.com/_/debian/
FROM debian:jessie

ENV JAVA_VERSION=8

# Install necessary packages including java 8 jre and clean up apt caches
RUN echo "deb http://ppa.launchpad.net/webupd8team/java/ubuntu trusty main" >> /etc/apt/sources.list.d/webupd8team-java.list && \
    echo "deb-src http://ppa.launchpad.net/webupd8team/java/ubuntu trusty main" >> /etc/apt/sources.list.d/webupd8team-java.list && \
    apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys EEA14886 && \
    echo debconf shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections && \
    echo debconf shared/accepted-oracle-license-v1-1 seen true | /usr/bin/debconf-set-selections

RUN apt-get update && apt-get --no-install-recommends install -y --force-yes \
    oracle-java${JAVA_VERSION}-installer && \ 
    apt-get clean autoclean && \
    apt-get autoremove -y && \
rm -rf /var/lib/{apt,dpkg,cache,log}/ /var/cache/oracle-jdk${JAVA_VERSION}-installer

RUN apt-get update && \
  apt-get install -y build-essential \
    wget \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    zip \
    python3-dev \
    python3-pip \
    r-base \
    pandoc \
    git \
    unzip

RUN mkdir -p /opt/software 

# Install fastQC:
RUN cd /opt/software && \
  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
  unzip fastqc_v0.11.8.zip && \
  chmod +x /opt/software/FastQC/fastqc
ENV PATH="/opt/software/FastQC:${PATH}"

# Install libStatGen and fastQValidator:
RUN cd /opt/software && \
  git clone https://github.com/statgen/libStatGen.git && \
  cd libStatGen && \
  make && \
  cd .. && \
  git clone https://github.com/statgen/fastQValidator.git && \
  cd fastQValidator && \
  make && \
make install

# Install samtools:
RUN cd /opt/software && \
  wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
  tar -xjf samtools-1.9.tar.bz2 && \
  cd samtools-1.9 && \
  ./configure --prefix=/opt/software/samtools && \
  make && \
  make install && \
  rm /opt/software/samtools-1.9.tar.bz2
ENV PATH="/opt/software/samtools/bin:${PATH}"

# Get the Picard JAR:
RUN mkdir -p /opt/software/picard && \
  cd /opt/software/picard && \
  wget https://github.com/broadinstitute/picard/releases/download/2.18.26/picard.jar
ENV PICARD_JAR=/opt/software/picard/picard.jar

# Install BWA
RUN mkdir -p /opt/software/ && \
  cd /opt/software/ && \
  wget https://github.com/lh3/bwa/archive/v0.7.17.tar.gz && \
  tar xzf v0.7.17.tar.gz && \
  cd bwa-0.7.17 && \
  make
ENV PATH="/opt/software/bwa-0.7.17:${PATH}"

# Get the GATK JAR:
RUN mkdir -p /opt/software/gatk && \
  cd /opt/software/gatk && \
  wget https://github.com/broadinstitute/gatk/releases/download/4.1.1.0/gatk-4.1.1.0.zip && \
  unzip gatk-4.1.1.0.zip
ENV GATK_JAR=/opt/software/gatk/gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar

# The scripts for performing pre-workflow checks:
RUN mkdir /opt/software/precheck
ADD check_fastq.py /opt/software/precheck/