FROM centos:centos5

RUN rpm -Uvh http://dl.fedoraproject.org/pub/epel/5/x86_64/epel-release-5-4.noarch.rpm
RUN yum install -y wget bzip2 bzip2-devel git gcc gcc-c++ patch zlib-devel make gcc44 gcc44-c++ cmake ncurses-devel unzip
RUN wget http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo
RUN wget --no-check-certificate https://copr.fedoraproject.org/coprs/praiskup/autotools/repo/epel-5/praiskup-autotools-epel-5.repo -O /etc/yum.repos.d/autotools.repo
RUN yum install -y devtoolset-2-gcc devtoolset-2-binutils devtoolset-2-gcc-c++ devtoolset-2-gcc-gfortran autotools-latest pkgconfig which

# install conda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -b -p /anaconda
ENV PATH=/opt/rh/devtoolset-2/root/usr/bin:/opt/rh/autotools-latest/root/usr/bin:/anaconda/bin:$PATH
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8
RUN conda update conda
