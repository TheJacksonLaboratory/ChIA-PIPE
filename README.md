## ChIA-PIPE  README
### The Jackson Laboratory for Genomic Medicine
### Version 1.0, (November 26, 2018) 
#

## (a) Initialize outer working directory
```bash
mkdir -p testing_chia_pipe
cd testing_chia_pipe
```

## (b) Clone ChIA-PIPE
```bash
git clone https://github.com/TheJacksonLaboratory/ChIA-PIPE.git
```

## *or* (b.ii) Download ChIA-PIPE
```bash
wget git@github.com:TheJacksonLaboratory/chia_pipe.zip
```

## (c) Install ChIA-PIPE dependencies
```bash
cd ChIA-PIPE/
dep_dir="dep_dir"
bash local_install_chia_pipe_dependencies.sh -i ${dep_dir}
```

## (d) Download test data from Dropbox
```bash
https://www.dropbox.com/sh/kiunzdfj74dpamh/AAB3u4vVGHIUSZoiPSjEYuYva?dl=0
mkdir -p fastq
cp  LDK0004-ds_*.fastq.gz  fastq
```

## (e) Launch ChIA-PIPE
```bash
bash chia_pipe-master/0.chia_pipe_shell.sh -c chia_pipe-master/example_config_file.sh
```

## *or* (e.ii) Launch ChIA-PIPE on HPC
```bash
qsub -F "--conf chia_pipe-master/example_config_file.sh" chia_pipe-master/0.chia_pipe_hpc.pbs
```
