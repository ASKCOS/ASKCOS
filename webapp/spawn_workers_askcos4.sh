source activate askcos
rm celery_logs/*.log
rm /home/mpiuser/askcosblade*

# Purge reservable queues to avoid weird conflicts
celery amqp queue.purge tb_c_worker
celery amqp queue.purge tb_c_worker_reservable
celery amqp queue.purge tb_coordinator
celery amqp queue.purge tb_coordinator_mcts
celery amqp queue.purge te_coordinator
celery amqp queue.purge sc_coordinator
celery amqp queue.purge ft_worker
celery amqp queue.purge tffp_worker
celery amqp queue.purge cr_coordinator
celery amqp queue.purge cr_network_worker

ssh mpiuser@18.4.94.201 "pkill -f 'celery'" || true 
ssh mpiuser@18.4.94.201 "cd /home/mpiuser/ASKCOS/ASKCOS_Website/; nohup spawn_scripts/spawn_scorers.sh > ~/askcosblade0$i.out 2> ~/askcosblade0$i.err < /dev/null"

for i in `seq 1 2`;
do
	ssh mpiuser@18.4.94.20$i "pkill -f 'celery'" || true 
	ssh mpiuser@18.4.94.20$i "cd /home/mpiuser/ASKCOS/ASKCOS_Website/; nohup spawn_scripts/spawn_scorers.sh > ~/askcosblade0$i.out 2> ~/askcosblade0$i.err < /dev/null"
done  

for i in `seq 3 3`;
do
	ssh mpiuser@18.4.94.20$i "pkill -f 'celery'" || true 
	ssh mpiuser@18.4.94.20$i "cd /home/mpiuser/ASKCOS/ASKCOS_Website/; nohup spawn_scripts/spawn_expanders_oldtree.sh > ~/askcosblade0$i.out 2> ~/askcosblade0$i.err < /dev/null"
done  

for i in `seq 4 7`;
do
	ssh mpiuser@18.4.94.20$i "pkill -f 'celery'" || true 
	ssh mpiuser@18.4.94.20$i "cd /home/mpiuser/ASKCOS/ASKCOS_Website/; nohup spawn_scripts/spawn_expanders.sh > ~/askcosblade0$i.out 2> ~/askcosblade0$i.err < /dev/null"
done  



source deactivate
