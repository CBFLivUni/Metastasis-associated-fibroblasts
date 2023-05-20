# sudo docker rm $(sudo docker ps -aq)
# sudo docker volume prune


docker run -d --cpus 30 -m 125g -p ****:**** -e ROOT=TRUE -e PASSWORD=****** --name ******* -e USERID=$(id -u) \
	-v $(pwd):/home/rstudio/project \
	-v /home/cbf-shared/CBF_projects:/home/rstudio/data\
	rocker/tidyverse:4.0.2

