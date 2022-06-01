FROM opensuse/leap:latest
LABEL version="1.1"
LABEL description="Dockerfile to build and run open-mmpol library"
RUN zypper --non-interactive in gcc \
                                gcc-c++ \ 
                                gcc-fortran \ 
                                make \
                                cmake \
                                python \
                                lapack-devel \
                                liblapack3 \ 
                                hdf5 \
                                hdf5-devel \
                                zlib-devel
RUN pip install ford
