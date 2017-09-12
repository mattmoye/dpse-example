FROM debian:7.11

## install dependencies
RUN apt-get update && apt-get install -y \
  wget \
  build-essential \
  pkg-config \
  gfortran \
  python-sympy \
  libatlas-base-dev

RUN wget --no-check-certificate https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.8.tgz \
    && echo "99180c51b0cde326a7e53f4661a5742baa1507aa  Ipopt-3.12.8.tgz" | sha1sum -c - \
    && mkdir -p /usr/src/Ipopt \
    && tar zxvf Ipopt-3.12.8.tgz -C /usr/src/Ipopt --strip-components=1 \
    && rm Ipopt-3.12.8.tgz

## insert your download URL from HSL here, or a URL to your private copy
## comment this out if you don't need the HSL 2011 solvers
RUN wget --no-check-certificate https://gracula.psyc.virginia.edu/public/software/coinhsl-2014.01.10.tar.gz \
    && mkdir -p /usr/src/Ipopt/ThirdParty/HSL/coinhsl \
    && tar zxvf coinhsl-2014.01.10.tar.gz -C /usr/src/Ipopt/ThirdParty/HSL/coinhsl --strip-components=1 \
    && rm coinhsl-2014.01.10.tar.gz

WORKDIR /usr/src/Ipopt
RUN cd ThirdParty/Metis \
    && ./get.Metis \
    && cd ../Mumps \
    && ./get.Mumps \
    && cd ../.. \
    && ./configure --with-blas="-lblas -llapack" --with-lapack="-llapack" --prefix="/usr/local" \
    && make -j $(nproc) \
    && make install \
    && ldconfig
#    && rm -rf /usr/src/Ipopt

WORKDIR /app
ADD . /app

## build the code using the example equations.txt so that it's in the image
RUN python model/makecode.py \
    && make


CMD ["/bin/bash"]
