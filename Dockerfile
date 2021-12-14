FROM alpine:3.11.5
#FROM alpine:20190925

ENV CFLAGS="-fPIC -O3"

RUN apk add wget git xz bzip2-static musl m4 autoconf tar xz-dev bzip2-dev build-base libpthread-stubs libzip-dev gfortran \
	    openssl-libs-static openblas-static pcre-dev curl llvm-dev curl-static bash

RUN mkdir -p /usr/local/include && \
    git clone --depth 1 https://github.com/ebiggers/libdeflate.git && \
    cd libdeflate && make -j4 CFLAGS="-fPIC -O3" install && \
    cd .. && rm -rf libdeflate && \
    git clone https://github.com/cloudflare/zlib cloudflare-zlib && \
    cd cloudflare-zlib && ./configure && make install && \
    cd .. && \
    rm -rf cloudflare-zlib

RUN cd / && \
    git clone -b v1.2.6 git://github.com/nim-lang/nim nim && \
    cd nim &&  \
    sh ./build_all.sh && \
    rm -rf csources && \
    echo 'PATH=/nim/bin:$PATH' >> ~/.bashrc && \
    echo 'PATH=/nim/bin:$PATH' >> ~/.bash_profile && \
    echo 'PATH=/nim/bin:$PATH' >> /etc/environment 

RUN apk add cmake openssl-dev && \
	wget https://libzip.org/download/libzip-1.6.1.tar.gz && \
	tar xzvf libzip-1.6.1.tar.gz && \
	cd libzip-1.6.1 && \
	mkdir build && cd build && \
	cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=/usr/local/ ../ && \
	make -j4 CFLAGS="-fPIC -O3" install && \
	cd ../../ && rm -rf libzip-1.6.1*

## libRmath.a

RUN apk update
RUN \
    apk add curl-dev

RUN \
    git clone https://github.com/SurajGupta/r-source/ && \
    cd r-source && chmod 755 ./configure && \
    ./configure -with-readline=no --with-x=no --disable-s3 --with-libdeflate && \
    cd src/nmath/standalone && make && \
    make install 

RUN \
    cd / && cp /r-source/src/nmath/standalone/libRmath.a /usr/local/lib/ && \
    cp /r-source/src/nmath/standalone/libRmath.so /usr/local/lib/
 
ENV PATH=:/root/.nimble/bin:/nim/bin/:$PATH	

RUN \
    git clone https://github.com/samtools/htslib && \
    cd htslib && git checkout 1.11 && autoheader && autoconf && \
    ./configure --disable-s3 --disable-libcurl --with-libdeflate && \
    make -j4 CFLAGS="-fPIC -O3" install && \
    cd ../ && \
    git clone https://github.com/samtools/bcftools && \
    cd bcftools && git checkout 1.10.2 && autoheader && autoconf && \
    ./configure --disable-s3 --disable-libcurl --with-libdeflate && \
    make -j4 CFLAGS="-fPIC -O3" install && \
    cd ../ && rm -rf htslib bcftools

ENV HTSLIB=system
ENV PATH=$PATH:~/.cargo/bin/

RUN \
    git clone https://gitlab.svi.edu.au/biocellgen-public/sgcocaller.git && \
    cd sgcocaller && \
    nimble install .
#ADD . /src/

#RUN nimble install -y https://gitlab.svi.edu.au/biocellgen-public/sgcocaller.git

#RUN ls ~/.nimble/lib/
ENV LD_LIBRARY_PATH=:/root/.nimble/lib/:/usr/local/lib/:$LD_LIBRARY_PATH

RUN cp /root/.nimble/bin/sgcocaller /usr/bin/

RUN /usr/bin/sgcocaller