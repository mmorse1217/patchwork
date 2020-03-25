# defines base image containing dependencies
FROM ubuntu:18.04 as patchwork-deps
RUN apt-get update && apt-get install -y \
    build-essential \ 
    cmake \ 
    git \
    extra-cmake-modules \
    freeglut3-dev \
    libblas-dev \
    liblapack-dev \
    mpich \
    wget \
    zlib1g-dev &&\
    apt-get purge -y curl && \
    apt-get autoremove -y && \
    #apt-get purge -y ca-certificates &&\
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install blendsurf
RUN mkdir libs/ && \
    cd libs &&\
    git clone https://github.com/mmorse1217/blendsurf &&\
    mkdir -p blendsurf/build/ && \
    cd blendsurf/build/ && \
    cmake -DCMAKE_MODULE_PATH=/usr/share/cmake-3.10/Modules/ ..  && \
    make 

# Install p4est 
RUN cd libs && \
    wget http://p4est.github.io/release/p4est-1.1.tar.gz && \
    tar xvf p4est-1.1.tar.gz && \
    rm p4est-1.1.tar.gz && \
    cd /libs/p4est-1.1 && \
    ./configure CC=mpicc F77=mpif77 FC=mpif90 --prefix=`pwd` --enable-mpi && \
    make && \
    make install 

# copy some relevant files into /patchwork required for blendsurf to run
RUN mkdir /patchwork && \
    cp /libs/blendsurf/ccsubmatall.dat /patchwork && \
    cp /libs/blendsurf/bdsurf_U_ONE.dat /patchwork

# specific location of blensurf and p4est
ENV BLENDSURF_DIR=/libs/blendsurf P4EST_DIR=/libs/p4est-1.1

# defines CI build: checks that patchwork core and renderer still compile
#FROM patchwork-deps as patchwork-build

# copy source code from repo into container 
#COPY . /patchwork/

# copy entrypoint.sh into container and execute its contents on "docker run"
#COPY entrypoint.sh /entrypoint.sh
#ENTRYPOINT ["/entrypoint.sh"]

#FROM patchwork-deps as patchwork-dev
WORKDIR /patchwork
CMD ["/bin/bash"]
