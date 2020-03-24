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

RUN apt-get update && apt-get install -y wget

# Install p4est 
RUN cd libs && \
    wget http://p4est.github.io/release/p4est-1.1.tar.gz && \
    tar xvf p4est-1.1.tar.gz && \
    rm p4est-1.1.tar.gz && \
    cd /libs/p4est-1.1 && \
    ./configure CC=mpicc F77=mpif77 FC=mpif90 --prefix=`pwd` --enable-mpi && \
    make && \
    make install 

RUN mkdir /patchwork
RUN cp /libs/blendsurf/ccsubmatall.dat /patchwork
RUN cp /libs/blendsurf/bdsurf_U_ONE.dat /patchwork

# defines CI build: checks that patchwork core and renderer still compile
FROM patchwork-deps as patchwork-build
#COPY . /patchwork/
#RUN mkdir -p /patchwork/build/ 
#WORKDIR /patchwork/build/
#ENV BLENDSURF_DIR=/libs/blendsurf P4EST_DIR=/libs/p4est-1.1
#RUN cmake -DCMAKE_MODULE_PATH=/usr/share/cmake-3.10/Modules/ ..  && \
#    make
#WORKDIR /patchwork/
COPY . /patchwork/
COPY entrypoint.sh /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]


#RUN CMAKE_PREFIX_PATH=/libs/blendsurf:/libs/p4est cmake -DCMAKE_MODULE_PATH=/usr/share/cmake-3.10/Modules/ .. && \
#    make
#        -DCOMPILE_RENDERER=True ..  && make 
#
## to use the container
FROM patchwork-deps as patchwork-dev

CMD ["/bin/bash"]
