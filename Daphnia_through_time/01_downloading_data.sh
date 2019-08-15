############################################################################
### Getting data from BGI ftp servers
############################################################################

wget <ftp download link for file or directory> --ftp-user=<username> --ftp-password=<password>

# Username: F19FTSEUET0079
# Password: DAPrecM_1

# Report
wget ftp://cdts-wh.genomics.cn/F19FTSEUET0079_DAPrecM/report.tar.gz --ftp-user=F19FTSEUET0079 --ftp-password=DAPrecM_1

# Readme
wget ftp://cdts-wh.genomics.cn/F19FTSEUET0079_DAPrecM/readme.txt --ftp-user=F19FTSEUET0079 --ftp-password=DAPrecM_1

# md5
wget ftp://cdts-wh.genomics.cn/F19FTSEUET0079_DAPrecM/md5.check --ftp-user=F19FTSEUET0079 --ftp-password=DAPrecM_1
wget ftp://cdts-wh.genomics.cn/F19FTSEUET0079_DAPrecM/md5.txt --ftp-user=F19FTSEUET0079 --ftp-password=DAPrecM_1

############################################################################

#!/bin/bash

#PBS -N getting_data
#PBS -l walltime=57:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

cat sample_names_BGI.txt | while read folder
do
    wget ftp://cdts-wh.genomics.cn/F19FTSEUET0079_DAPrecM/Clean/${folder}/* \
    --ftp-user=F19FTSEUET0079 --ftp-password=DAPrecM_1
done




