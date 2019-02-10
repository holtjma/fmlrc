#!/bin/bash
#Download this script via "git clone https://github.com/holtjma/fmlrc.git"
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

#download short-read ecoli
if [ ! -f s_6_1.fastq.gz ]; then
    curl -o s_6_1.fastq.gz http://spades.bioinf.spbau.ru/spades_test_datasets/ecoli_mc/s_6_1.fastq.gz
fi
if [ ! -f s_6_2.fastq.gz ]; then
    curl -o s_6_2.fastq.gz http://spades.bioinf.spbau.ru/spades_test_datasets/ecoli_mc/s_6_2.fastq.gz
fi

#download long-read ecoli and convert to fasta format
if [ ! -f PacBioCLR/PacBio_10kb_CLR.fasta ]; then
    curl -o Ecoli_MG1655_pacBioToCA.tgz http://files.pacb.com/datasets/secondary-analysis/e-coli-k12-de-novo/1.3.0/Ecoli_MG1655_pacBioToCA.tgz
    tar -xvzf Ecoli_MG1655_pacBioToCA.tgz
    awk 'NR%4==1||NR%4==2' ./PacBioCLR/PacBio_10kb_CLR.fastq | tr "@" ">" > ./PacBioCLR/PacBio_10kb_CLR.fasta
fi

#build the bwt
if [ ! -f ./ecoli_comp_msbwt.npy ]; then
    mkdir temp
    gunzip -c s_6_?.fastq.gz | awk "NR % 4 == 2" | sort -T ./temp | tr NT TN | ./ropebwt2/ropebwt2 -LR | tr NT TN | ../fmlrc-convert ./ecoli_comp_msbwt.npy
fi

#run fmlrc
NUM_PROCS=4
../fmlrc -p $NUM_PROCS -e 400 ./ecoli_comp_msbwt.npy ./PacBioCLR/PacBio_10kb_CLR.fasta ./corrected_final.fa