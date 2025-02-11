FROM rocker/r-ver:4.0.3

RUN apt-get update --fix-missing -qq && \
	apt-get install -y -q \
	vim \
	git \
	time \
	python3 \
	python3-pip \
	libz-dev \
	libcurl4-gnutls-dev \
	libxml2-dev \
	libssl-dev \
	libpng-dev \
	libjpeg-dev \
	libbz2-dev \
	liblzma-dev \
	libncurses5-dev \
	libncursesw5-dev \
	libgl-dev \
	libgsl-dev \
	sysstat \
	watch \
	&& apt-get clean \
	&& apt-get purge \
	&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN cd /tmp && git clone --branch 1.11.0 https://github.com/samtools/htslib.git && cd htslib && make && make install
RUN cd /tmp && git clone --branch 1.11 https://github.com/samtools/samtools.git && cd samtools && make && make install
RUN cd /tmp && git clone --branch 1.11 https://github.com/samtools/bcftools.git && cd bcftools && make && make install

RUN pip3 install numpy cython pandas
RUN pip3 install pysam
RUN pip3 install HTSeq
RUN pip3 install scikit-learn statsmodels

RUN R -e 'install.packages(c("BiocManager", "devtools", "argparse", "dbplyr"))'
RUN R -e 'BiocManager::install("tximport")'
RUN R -e 'BiocManager::install("rhdf5")'
RUN R -e 'BiocManager::install("GenomicFeatures")'
RUN R -e 'BiocManager::install("EnrichmentBrowser")'
RUN R -e 'BiocManager::install("BANDITS")'
RUN R -e 'BiocManager::install("DRIMSeq")'
RUN R -e 'BiocManager::install("DEXSeq")'


ADD data /home/data
ADD scripts /home/scripts
ADD software /home/software

RUN cd /home/software && git clone --branch v2.3 https://github.com/comprna/SUPPA.git
