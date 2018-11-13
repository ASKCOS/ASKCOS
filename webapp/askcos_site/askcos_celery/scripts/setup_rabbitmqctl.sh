sudo rabbitmqctl add_user ccoley password
sudo rabbitmqctl add_vhost askcos_vhost
sudo rabbitmqctl set_user_tags ccoley ccoley
sudo rabbitmqctl set_permissions -p askcos_vhost ccoley ".*" ".*" ".*"
