FROM daschapps:latest AS builder

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
   bison \
   build-essential \
   ca-certificates \
   curl \
   giza-dev \
   libatlas-base-dev \
   libbz2-dev \
   libcurl4-openssl-dev \
   libfftw3-dev \
   libgd-dev \
   libmariadb-dev \
   libopenmpi-dev \
   libplplot-dev \
   libtool \
   octave-dev \
   pkg-config \
   python3 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# daschapps base image provides wcstools/libwcs

COPY install_funtools.sh /tmp
RUN /tmp/install_funtools.sh && rm -f /tmp/install_funtools.sh

COPY install_starbase.sh starbase-fixes.patch /tmp
RUN /tmp/install_starbase.sh && rm -f /tmp/install_starbase.sh /tmp/starbase-fixes.patch

COPY install_sextractor.sh sextractor-fixes.patch /tmp
RUN /tmp/install_sextractor.sh && rm -f /tmp/install_sextractor.sh /tmp/sextractor-fixes.patch

COPY install_scamp.sh scamp-fixes.patch /tmp
RUN /tmp/install_scamp.sh && rm -f /tmp/install_scamp.sh /tmp/scamp-fixes.patch

COPY install_astrometry_net.sh /tmp
RUN /tmp/install_astrometry_net.sh && rm -f /tmp/install_astrometry_net.sh

RUN mkdir /build
COPY Makefile /build
COPY bin_scripts /build/bin_scripts
COPY scripts /build/scripts
COPY src /build/src

RUN cd /build \
 && mkdir bin \
 && make clean all install PREFIX=/dasch

FROM daschapps:latest

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
   ca-certificates \
   emacs-nox \
   gdb \
   imagemagick \
   iputils-ping \
   libatlas3-base \
   libbz2-1.0 \
   libcpgplot0 \
   libcurl4 \
   libfftw3-double3 \
   libgd3 \
   libmariadb3 \
   libopenmpi3 \
   libplplot17 \
   mysql-client \
   octave \
   octave-image \
   octave-io \
   octave-optim \
   octave-statistics \
   octave-struct \
   openmpi-bin \
   plplot-driver-cairo \
   python-is-python3 \
   python3 \
   python3-astropy \
   python3-matplotlib \
   strace \
   telnet \
   tcsh \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /dasch /dasch

COPY etc /dasch/etc

ENV DASCH_SCRIPTS=/dasch/pipeline
