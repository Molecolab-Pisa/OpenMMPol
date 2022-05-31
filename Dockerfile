FROM opensuse/leap:latest
LABEL version="1.0"
LABEL description="Dockerfile to build and run open-mmpol library"

# Preparing environment for installing PGI compiler
RUN zypper --non-interactive in gcc \
                                gcc-c++ \ 
                                gcc-fortran \ 
                                make \
                                python \
                                lapack-devel \
                                liblapack3
