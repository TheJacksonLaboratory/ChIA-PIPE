###############################################################
# Draft Dockerfile to build ChIA-PET Utilities (CPU)
# ChIA-PET Utilities provides modules with core data-processing 
# functionality for ChIA-PET data
# Based on "CentOS release 6.5 (Final)"
###############################################################

# Set the base image to CentOS6
FROM centos:6.6

# File Author / Maintainer
MAINTAINER Daniel Capurso

## Dependencies
# zlib v 1.2.8
# ChiaSigScaled v 19.4.17n
# git
# bwa       (Compiled during CPU make ??)
# BWT-SW    (Compiled during CPU make ??)


# Update the repository sources list
# Install base packages: git, python, wget, unzip, R
RUN yum -y update && yum -y install \
    git \
    curl \
    tar \
    unzip \
    make \
    automake \
    gcc \
    gcc-c++ \
    kernel-devel

## Make outer working directory
RUN mkdir cpu-dir

## Download and compile zlib (v 1.2.8)
RUN cd cpu-dir && curl -O https://www.zlib.net/fossils/zlib-1.2.8.tar.gz
RUN cd cpu-dir && tar -xzvf zlib-1.2.8.tar.gz
RUN cd cpu-dir/zlib-1.2.8 && ./configure && make 

## Download and compile Chiasig
RUN cd cpu-dir && curl -O http://folk.uio.no/jonaspau/chiasig/ChiaSigCPPv093.zip
RUN cd cpu-dir && unzip ChiaSigCPPv093.zip
RUN cd cpu-dir/ChiaSigCPPv093 && make

## Compile CPU
#RUN cd cpu-dir && make
