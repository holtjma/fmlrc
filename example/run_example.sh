#!/bin/bash
#Download this repository via "git clone https://github.com/holtjma/fmlrc.git"
#cd to "fmlrc/examples" and run "./run_example.sh"

#download and build ropebwt2
if [ ! -f ./ropebwt2/ropebwt2 ]; then
    git clone https://github.com/lh3/ropebwt2.git
    cd ropebwt2; make; cd ..
fi

#download and build fmlrc
if [ ! -f ../fmlrc ]; then
    cd ..; make; cd example
fi

DATADIR="example1"
if [ ! -d ${DATADIR} ]; then
    mkdir ${DATADIR}
fi

#download short-read ecoli
if [ ! -f ${DATADIR}/ERR022075_1.fastq.gz ]; then
    curl -o ${DATADIR}/ERR022075_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022075/ERR022075_1.fastq.gz
fi
if [ ! -f ${DATADIR}/ERR022075_2.fastq.gz ]; then
    curl -o ${DATADIR}/ERR022075_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022075/ERR022075_2.fastq.gz
fi

#download long-read ecoli and convert to fasta format
if [ ! -f ${DATADIR}/PacBioCLR/PacBio_10kb_CLR.fasta ]; then
    curl -o ${DATADIR}/Ecoli_MG1655_pacBioToCA.tgz http://files.pacb.com/datasets/secondary-analysis/e-coli-k12-de-novo/1.3.0/Ecoli_MG1655_pacBioToCA.tgz
    tar -xvzf ${DATADIR}/Ecoli_MG1655_pacBioToCA.tgz -C ${DATADIR}
    awk 'NR%4==1||NR%4==2' ${DATADIR}/PacBioCLR/PacBio_10kb_CLR.fastq | tr "@" ">" > ${DATADIR}/PacBioCLR/PacBio_10kb_CLR.fasta
fi

#build the bwt
if [ ! -f ${DATADIR}/ecoli_comp_msbwt.npy ]; then
    mkdir temp
    gunzip -c ${DATADIR}/ERR022075_?.fastq.gz | awk "NR % 4 == 2" | sort -T ./temp | tr NT TN | ./ropebwt2/ropebwt2 -LR | tr NT TN | ../fmlrc-convert ${DATADIR}/ecoli_comp_msbwt.npy
fi

#run fmlrc
NUM_PROCS=4
../fmlrc -p $NUM_PROCS -e 400 ${DATADIR}/ecoli_comp_msbwt.npy ${DATADIR}/PacBioCLR/PacBio_10kb_CLR.fasta ${DATADIR}/corrected_final.fa