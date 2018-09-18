# Set the directory for local installation of ChIA-PIPE dependencies
# This directory should be in your $PATH
install_dir="/home/capurd/bin"

## Move to the install directory
cd ${install_dir}

## Install pigz
wget http://zlib.net/pigz/pigz-2.4.tar.gz
tar -xzvf pigz-2.4.tar.gz
cd pigz-2.4
make
cp pigz unpigz ../
cd ../

## Install java/1.7.0


## Install perl/5.26.0


## Install bedtools/2.26.0


## Install samtools/1.5


## Install R/3.2.1


## Install MACS/2.1.0.20151222






