apt-get install dos2unix #not very necessary

#set up a conda at your own space first, I am not going to deal with the sudo
#follow https://docs.conda.io/projects/conda/en/stable/user-guide/install/linux.html
source /home/rxr456/miniconda3/etc/profile.d/conda.sh
#conda automatically installs a version of Python alongside it

#I like salmon!
apt-get install salmon #or module load salmon
sudo salmon index -t /home/rxr456/Homo_sapiens.GRCh38.cdna.all.fa.gz -i /home/rxr456/salmon_index/transcripts_index -k 31 -p 48

#RSEM
module load R
module load STAR

#RSEM relies on perl, get perl locally
curl -L http://xrl.us/installperlnix | bash
echo 'export PATH=/home/rxr456/perl5/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
echo 'export PERL5LIB=/home/rxr456/perl5/lib/perl5:$PERL5LIB' >> ~/.bashrc
source ~/.bashrc

wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
gunzip Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
gunzip Homo_sapiens.GRCh38.112.gtf.gz

cd RSEM
wget https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz
tar -xzf https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz
cd RSEM-1.3.3
make
export PATH=$PATH:/home/rxr456/RSEM/RSEM-1.3.3

#run stuff below one more time incase something missing
module load Java
module load R
module load STAR
module load GCC
export PATH=$PATH:/home/rxr456/RSEM/RSEM-1.3.3 
echo 'export PATH=/home/rxr456/perl5/bin:$PATH' >> ~/.bashrc
source ~/.bashrc  
echo 'export PERL5LIB=/home/rxr456/perl5/lib/perl5:$PERL5LIB' >> ~/.bashrc
source ~/.bashrc
source /home/rxr456/miniconda3/etc/profile.d/conda.sh #make sure you have the conda and python correctly installed

#Using rsem-prepare-reference, you prepare the reference data for both RSEM and STAR.
rsem-prepare-reference --gtf /home/rxr456/Homo_sapiens.GRCh38.112.gtf /home/rxr456/Homo_sapiens.GRCh38.dna.primary_assembly.fa /home/rxr456/RSEM_ref/RSEM_ref --star