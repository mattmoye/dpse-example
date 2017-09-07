FROM debian:7.11

## install dependencies
RUN apt-get update && apt-get install -y \
  wget \
  build-essential \
  gfortran \
  python-pip \
  python-sympy \
  libatlas-base-dev

RUN wget --no-check-certificate https://www.coin-or.org/download/source/Ipopt/Ipopt-3.10.2.tgz \
    && echo "6a19cac8772c050bd52214a280b69b96ad18b6fa  Ipopt-3.10.2.tgz" | sha1sum -c - \
    && mkdir -p /usr/src/Ipopt \
    && tar zxvf Ipopt-3.10.2.tgz -C /usr/src/Ipopt --strip-components=1 \
    && rm Ipopt-3.10.2.tgz

## insert your download URL from HSL here, or a URL to your private copy
## comment this out if you don't need the HSL 2011 solvers
RUN wget --no-check-certificate https://gracula.psyc.virginia.edu/public/software/coinhsl-2014.01.10.tar.gz \
    && mkdir -p /usr/src/coinhsl \
    && tar zxvf coinhsl-2014.01.10.tar.gz -C /usr/src/coinhsl --strip-components=1 \
    && find /usr/src/coinhsl/ -name "m*.f" -exec cp {} /usr/src/Ipopt/HSL ';' \
    && rm coinhsl-2014.01.10.tar.gz

WORKDIR /usr/src/Ipopt
RUN cd ThirdParty/Metis \
    && ./get.Metis \
    && cd ../Mumps \
    && ./get.Mumps \
    && cd ../.. \
    && ./configure --with-blas="-lblas -llapack" --with-lapack="-llapack" --prefix="/usr/local" \
    && make \
    && make install
    && apt-get purge -y --auto-remove build-essential \
    && rm -rf /usr/src/Ipopt /usr/src/coinshsl

WORKDIR /app
ADD . /app

## build the code using the example equations.txt so that it's in the image
RUN python models/makecode.py \
    && make \
    && biohh1


CMD ["/bin/bash"]
