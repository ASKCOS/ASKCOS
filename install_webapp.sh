### Meant for Python 2.7.6 on Ubuntu 16.04 running on an AWS EC2 instance

sudo add-apt-repository -y ppa:chris-lea/redis-server 
sudo apt-get update
sudo apt-get install redis-server
sudo service redis start 

echo "deb https://dl.bintray.com/rabbitmq/debian xenial main" | sudo tee /etc/apt/sources.list.d/bintray.rabbitmq.list 
wget -O- https://dl.bintray.com/rabbitmq/Keys/rabbitmq-release-signing-key.asc |
     sudo apt-key add - 
sudo apt-get update 


wget https://packages.erlang-solutions.com/erlang-solutions_1.0_all.deb 
sudo dpkg -i erlang-solutions_1.0_all.deb
sudo apt-get update
sudo apt-get install erlang 

sudo apt-get install rabbitmq-server # on OS X, brew install rabbitmq

sudo service redis start 
sudo service rabbitmq-server start

# change ASKCOS_Website settings for celery.py

pip install uwsgi

udo apt-get install nginx
sudo /etc/init.d/nginx start 

# test that it is working with uwsgi --http :8000 --wsgi-file wsgi.py
# test by curl http://localhost:8000

# In website folder
wget https://raw.githubusercontent.com/nginx/nginx/master/conf/uwsgi_params

# Edit /etc/nginx/nginx.conf according to http://uwsgi-docs.readthedocs.io/en/latest/tutorials/Django_and_nginx.html

#Start the server

# Run wsgi with uwsgi --socket :8000 --wsgi-file wsgi.py
sudo service nginx start

# Test from outside (might need to open up security settings)

# Test Celery part of website by loading context rec workers:
celery -A askcos_site worker -c 1 -Q cr_coordinator -n "cr_coordinator@$(hostname)" --max-tasks-per-child 1000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 1 -Q cr_network_worker -n "cr_network_worker@$(hostname)" --max-tasks-per-child 5000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
# On the website, go to /context/ and make a request for the Neural Network context recommender
# Note: these workers should only take ~ 500 MB of RAM, so they can be tested on the main EC2 instance

