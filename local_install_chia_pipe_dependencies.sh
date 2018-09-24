## Perform local installation of dependencies for ChIA-PIPE

## The help message:
function usage
{
    echo -e "usage: bash local_install_chia_pipe_dependencies.sh -i INSTALL_DIR
    " 
}

# You will need to have write permission to the local installation directory.
# The local installation directory should be in your $PATH.

## Example value for the argument:
# install_dir="/path/for/installation/"


## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -i | --install_dir )    shift
                                install_dir=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done


## Move to the install directory
cd ${install_dir}

## Insert install directory at front of PATH
export PATH="${install_dir}:${PATH}"


## Install pigz
wget http://zlib.net/pigz/pigz-2.4.tar.gz
tar -xzvf pigz-2.4.tar.gz
cd pigz-2.4
make
cp pigz unpigz ../
cd ../
rm -r pigz-2.4


## Install java/1.8
wget -c --header "Cookie: oraclelicense=accept-securebackup-cookie" \
http://download.oracle.com/otn-pub/java/jdk/\
8u131-b11/d54c1d3a095b4ff2b6607d096fa80163/jdk-8u131-linux-x64.tar.gz
#
tar -xzvf jdk-8u131-linux-x64.tar.gz
cp -r jdk1.8.0_131/jre/ .
yes | rm -r jdk1.8.0_131 
ln -s jre/bin/java java


## Install perl/5.26.0
wget http://www.cpan.org/src/5.0/perl-5.26.0.tar.gz
tar -xzvf perl-5.26.0.tar.gz
cd perl-5.26.0
./Configure -des -Dprefix=${install_dir}
make
make test
make install
cd ..
yes | rm -r perl-5.26.0


## Install bedtools/2.26.0
wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/\
bedtools-2.26.0.tar.gz
#
tar -xzvf bedtools-2.26.0.tar.gz
cd bedtools2
make
cd bin
cp * ../../
cd ../../
rm -r bedtools2


## Install conda
wget https://repo.anaconda.com/archive/Anaconda2-5.2.0-Linux-x86_64.sh
bash Anaconda2-5.2.0-Linux-x86_64.sh
# ... Answer questions and set install dir to:
# /path/to/install/dir/anaconda2
#
ln -s anaconda2/bin/python python
ln -s anaconda2/bin/conda conda


## Install pysam (Python package, not included with conda)
./conda install -c bioconda pysam


## Install MACS peak caller
./conda install -c bioconda macs2
ln -s anaconda2/bin/macs2 macs2


## Install samtools/1.5
wget https://sourceforge.net/projects/samtools/files/samtools/1.5/\
samtools-1.5.tar.bz2
#
tar -xvjf samtools-1.5.tar.bz2
cd samtools-1.5
./configure --disable-lzma
make
cp samtools ../
cd ../
rm -r samtools-1.5


## Install R/3.2.1
http://lib.stat.cmu.edu/R/CRAN/src/base/R-3/R-3.2.1.tar.gz
tar -xzvf R-3.2.1.tar.gz
cd R-3.2.1
./configure --prefix=${install_dir}
make
cd ../
ln -s R-3.2.1/bin/R R
