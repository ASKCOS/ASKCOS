sudo service apache2 stop
source activate askcos
conda uninstall cairo -y
conda install cairo -y
conda install matplotlib -y
sudo service apache2 start
