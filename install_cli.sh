# (1) Set up miniconda environment
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p $HOME/miniconda
export PATH=~/miniconda/bin:$PATH
echo 'export PATH=~/miniconda/bin:$PATH' >> ~/.bashrc

# (2) Prepare system dependencies before setting up python packages
sudo apt install heimdal-dev
sudo apt install libkrb5-dev 
sudo apt-get install build-essential

# (3) Create "askcos" conda environment
conda-env create -f askcos.yml -n askcos

# (4) Obtain code [NOTE: CANNOT CLONE PRIVATE REPOS WITHOUT ENTERING USERNAME/PASSWORD]
mkdir ASKCOS
cd ASKCOS
git clone https://github.com/connorcoley/RDChiral
git clone https://github.com/connorcoley/SCScore
#git clone https://github.com/connorcoley/Make-It # private repo
#git clone https://github.com/connorcoley/ASKCOS_Website # private repo
export PYTHONPATH=~/ASKCOS/RDChiral:~/ASKCOS/SCScore:~/ASKCOS/Make-It:~/ASKCOS/ASKCOS_Website:$PYTHONPATH
echo 'export PYTHONPATH=~/ASKCOS/RDChiral:~/ASKCOS/SCScore:~/ASKCOS/Make-It:~/ASKCOS/ASKCOS_Website:$PYTHONPATH' >> ~/.bashrc 

# (5) [OPTIONAL] Create a link between Make-It/makeit/data folder and wherever it is stored
cd ..
rm -r ASKCOS/Make-It/makeit/data
ln -s /mnt/data ASKCOS/Make-It/makeit/data

# (6) Change data path names in Make-It/makeit/global_config.py and ASKCOS_Website/askcos_site/settings.py 
# note: there are some hard links to the Mongo database currently

# (7) Install RDKit and pymongo
source activate askcos
conda install rdkit 
pip install pymongo
sudo apt install libxrender-dev 

export PYTHONPATH=~/miniconda/envs/askcos/lib/python2.7/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=~/miniconda/envs/askcos/lib:$LD_LIBRARY_PATH

echo 'export PYTHONPATH=~/miniconda/envs/askcos/lib/python2.7/site-packages:$PYTHONPATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=~/miniconda/envs/askcos/lib:$LD_LIBRARY_PATH' >> ~/.bashrc

# (8) Install more recent version of tensorflow
pip install tensorflow==1.4.1

# (9) Set keras backend default to theano by editing ~/.keras/keras.json

# (10) Test make-it running (note: still req connection to askcos2 mongo db)
python ASKCOS/Make-It/makeit/application/run.py --TARGET atropine
