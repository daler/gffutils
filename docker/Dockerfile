FROM ubuntu:14.04
RUN apt-get update && apt-get install -y wget build-essential
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -b -p /anaconda
RUN apt-get -y install git
ENV PATH=/anaconda/bin:$PATH
RUN conda update conda
