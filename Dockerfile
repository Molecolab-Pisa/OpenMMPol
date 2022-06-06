FROM opensuse/leap:latest
LABEL version="1.1"
LABEL description="Dockerfile to build and run open-mmpol library"
RUN zypper --non-interactive in gcc \
                                gcc-c++ \ 
                                gcc-fortran \ 
                                make \
                                cmake \
                                python3-pip \
                                python3 \
                                lapack-devel \
                                liblapack3 \ 
                                hdf5 \
                                hdf5-devel \
                                zlib-devel
RUN pip install ford
RUN sed -i -e 's/subprocess.run(command, check=True, capture_output=True, text=True)/subprocess.run(command, check=True)/g' /usr/lib/python3.6/site-packages/ford/__init__.py 
