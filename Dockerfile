FROM opensuse/leap:15.6
LABEL version="1.8"
LABEL description="Dockerfile to build and run open-mmpol library"
RUN zypper --non-interactive install \
                                cJSON-devel \
                                cmake \
                                curl \
                                gcc \
                                gcc-c++ \
                                gcc-fortran \
                                gcovr \
                                git \
                                gzip \
                                hdf5 \
                                hdf5-devel \
                                kernel-devel \
                                kmod \
                                lapack-devel \
                                lcov \
                                liblapack3 \
                                make \
                                procps \
                                python311-numpy \
                                python311-pip \
                                python311-pybind11 \
                                python311-pybind11-devel \
                                python311-virtualenv \
                                python311 \
                                tar \
                                wget \
                                zlib-devel
# lcov_cobertura needed to convert lcov output to cobertura format (needed by gitlab)
RUN pip install lcov_cobertura
# Installation of Ford (needed for documentation)
RUN pip install -Iv ford==6.1.11
#Intel Compilers suite
# RUN zypper --non-interactive addrepo https://yum.repos.intel.com/oneapi oneAPI
# RUN curl https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB > intel.key
# RUN rpm --import intel.key
# RUN zypper --non-interactive --gpg-auto-import-keys install intel-basekit \
#                                                             intel-hpckit
# RUN wget https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-1_12_2.tar.gz; \
#     tar xvf hdf5-1_12_2.tar.gz; \
#     rm hdf5-1_12_2.tar.gz; \
#     cd hdf5-hdf5-1_12_2; \
#     source /opt/intel/oneapi/setvars.sh; \
#     CC=icx CXX=icpx FC=ifx ./configure --prefix /opt/intel/hdf5-1.12.2 --enable-fortran --enable-build-mode=production --enable-shared; \
#     make; make install; \
#     cd -; rm -rf hdf5-1_12_2.tar.gz;
# NVCompilers suite
# RUN zypper --non-interactive addrepo https://developer.download.nvidia.com/hpc-sdk/sles/nvhpc.repo
# RUN zypper --non-interactive --gpg-auto-import-keys --no-gpg-checks install nvhpc
# RUN wget https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-1_12_2.tar.gz; \
#     tar xvf hdf5-1_12_2.tar.gz; \
#     rm hdf5-1_12_2.tar.gz; \
#     cd hdf5-hdf5-1_12_2; \
#     export PATH=/opt/nvidia/hpc_sdk/`uname -s`_`uname -m`/2024/compilers/bin:$PATH; \
#     CC=nvc CXX=nvcc FC=nvfortran ./configure --prefix /opt/nvidia/hdf5-1.12.2 --enable-fortran --enable-build-mode=production --enable-shared; \
#     make; make install; \
#     cd -; rm -rf hdf5-1_12_2.tar.gz;
