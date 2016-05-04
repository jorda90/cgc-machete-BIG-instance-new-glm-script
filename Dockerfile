FROM centos:centos6
#comes with python/2.7.5
MAINTAINER Gillian Lee Hsieh <glhsieh@stanford.edu>
#RUN yum update -y && yum install -y wget git gcc unzip gcc-c++ zlib-devel openssl-devel sqlite-devel bzip2-devel ncurses-devel lapack-dev blas-dev
RUN yum update -y && yum groupinstall -y 'Development Tools' && yum install -y wget zlib-devel openssl-devel sqlite-devel bzip2-devel ncurses-devel lapack-dev blas-dev
#lapack-dev blas-dev for installing scipy in Python
#Development Tools installs 28 packages, and their dependencies. The number of dependencies installed on a base image of centos:centos6 were 101, notably of which are Perl v5.10.1, git v1.7.1, unzip
ENV DATA=/home/data IndelIndices=IndelIndices HG19exons=HG19exons circularRNApipeline_Standalone=circularRNApipeline_Standalone
RUN mkdir /src 
#gcc-c++ needed for running g++ to make tools, such as Bowtie2
#gcc needed for installing Python
RUN git clone https://github.com/nathankw/KNIFE.git /src/knife && \
	cd /src/knife && \	
	git remote add upstream https://github.com/lindaszabo/KNIFE.git && \
	git pull
RUN git clone https://github.com/nathankw/MACHETE.git /src/machete && \
	cd /src/machete && \
	git remote add upstream https://github.com/gillianhsieh/MACHETE && \
	git pull
#INSTALL Python 2.7.10
RUN mkdir -p /src/Python /src/software/Python && \
	cd /src/Python && \
	wget https://www.python.org/ftp/python/2.7.10/Python-2.7.10.tgz && \
	tar -zxf Python-2.7.10.tgz && \
	cd Python-2.7.10 && \
	./configure --prefix=/src/software/Python && \
	make && \
	make install && \
	export PATH=/src/software/Python/bin:${PATH} && \
	wget https://bootstrap.pypa.io/get-pip.py && \
	python get-pip.py && \
	pip install scipy
#INSTALL TBB (Threading Building Blocks) from Intel
#Needed for intalling Bowtie1 and Bowtie2 with parallelism enabled (to use the -p argument).
RUN mkdir /src/TBB && \
	cd /src/TBB && \
	wget https://www.threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb44_20160128oss_src_0.tgz && \
	tar -zxf tbb44_20160128oss_src_0.tgz && \
	cd tbb44_20160128oss && \
	gmake && \
	. build/linux_intel64_gcc_cc4.4.7_libc2.12_kernel4.1.19_release/tbbvars.sh
#INSTALL Bowtie1.1.1
RUN mkdir /src/Bowtie1 && \
	cd /src/Bowtie1 && \
	wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-src.zip && \
	unzip bowtie-1.1.2-src.zip && \
	cd bowtie-1.1.2 && \
#	make WITH_TBB=1 && \
	make && \
	make install
	
#INSTALL Bowtie2 2.2.8
RUN mkdir /src/Bowtie2 && \
	cd /src/Bowtie2 && \
	wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.8/bowtie2-2.2.8-source.zip && \
	unzip bowtie2-2.2.8-source.zip && \
	cd bowtie2-2.2.8 && \
#	make WITH_TBB=1 && \
	make && \
	make install
#INSTALL R >= 3.0.2
#As stated in https://cran.r-project.org/bin/linux/redhat/README: 
#	The Fedora RPMs for R have been ported to RHEL by the project Extra Packages for Enterprise Linux (EPEL)
RUN rpm -Uvh http://download.fedoraproject.org/pub/epel/6/i386/epel-release-6-8.noarch.rpm && yum install -y R
#Installs R v3.2.3

#INSTALL samtools/1.3. Needed for knife.
RUN mkdir /src/samtools && \
	cd /src/samtools && \
	wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 && \
	tar -jxf samtools-1.3.tar.bz2 && \
	cd samtools-1.3 && \
	./configure --without-curses && \
	make && \
	make install	

#perl installation not neccessary since Development Tools, which yum installed earlier, includes Perl v5.10.1.
#INSTALL Perl. Needed for knife. On sherlock v5.10.1, I'll grab the latest, however.
#RUN mkdir /src/Perl && \
#		cd /src/Perl && \
#		wget http://www.cpan.org/src/5.0/perl-5.22.1.tar.gz && \
#		tar -zxf perl-5.22.1.tar.gz && \
#		cd perl-5.22.1 && \
#		sh Configure -de && \
#		make && 
#		make install &&

##### ADD KNIFE Data Dependencies
#Note: The directory itself is not copied, just its contents.
#If <dest> does not end with a trailing slash, it will be considered a regular file and the contents of <src> will be written at <dest>.
#If <dest> doesn’t exist, it is created along with all missing directories in its path.
ADD /${circularRNApipeline_Standalone} ${DATA}/${circularRNApipeline_Standalone}/

#### ADD MACHETE Data Dependencies
#ADD HG19exons. Location of HG19exons was formerly called PICKLEDIR
ADD /${HG19exons} ${DATA}/${HG19exons}/

#ADD REG_INDEL_INDICES
RUN mkdir ${DATA}/${IndelIndices}
ADD /${IndelIndices} ${DATA}/${IndelIndices}/

ENTRYPOINT []
LABEL version="1.0" description="Detects gene fusions"
