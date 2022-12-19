FROM opensuse/leap:latest
LABEL version="1.2"
LABEL description="Dockerfile to build and run open-mmpol library"
RUN zypper --non-interactive in gcc \
                                gcc-c++ \ 
                                gcc-fortran \ 
                                make \
                                cmake \
                                python-pybind11-common-devel \
                                python3-numpy \
                                python3-pip \
                                python3-pybind11 \
                                python3-pybind11-devel \
                                python3 \
                                lapack-devel \
                                liblapack3 \ 
                                hdf5 \
                                hdf5-devel \
                                valgrind \
                                zlib-devel
RUN pip install -Iv ford==6.1.11
RUN sed -i -e 's/subprocess.run(command, check=True, capture_output=True, text=True)/subprocess.run(command, check=True)/g' /usr/lib/python3.6/site-packages/ford/__init__.py
RUN zypper --non-interactive addrepo https://yum.repos.intel.com/oneapi oneAPI
RUN zypper --non-interactive --gpg-auto-import-keys install intel-basekit \
                                                            intel-hpckit
RUN zypper --non-interactive addrepo https://developer.download.nvidia.com/hpc-sdk/sles/nvhpc.repo
RUN zypper --non-interactive --gpg-auto-import-keys install nvhpc
