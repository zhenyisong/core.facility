#---
# @author  Yisong Zhen
# @since   2018-01-31
# @update  2018-01-31
# Undemy course, thanks.
# install docker
# install edge version
#----


#---
# docker need root power, otherwise, cannot install or use it.
# I tried using otherway.
# https://docs.docker.com/engine/installation/linux/docker-ce/centos/#os-requirements
# first we should check the system information
# on AWS. I tried it, it seems to be centos7
# https://github.com/docker/for-linux/issues/20
# Docker for Centos - Install Instructions Causing Conflicts
#---
cat /etc/os-release

sudo yum remove docker \
                docker-common \
                docker-selinux \
                docker-engine

sudo yum install -y yum-utils \
  device-mapper-persistent-data \
  lvm2

sudo yum-config-manager \
    --add-repo \
    https://download.docker.com/linux/centos/docker-ce.repo


sudo yum-config-manager --enable docker-ce-edge
sudo yum makecache fast
sudo yum -y --enablerepo=rhui-REGION-rhel-server-extras install container-selinux

sudo yum install docker-ce


# install docker 
curl -sSL https://get.docker.com/ | sh


# install docker machine
# install docker docker compose
#
sudo usermod -aG docker zhenyisong
sudo docker version

sudo service docker stop
sudo service docker start

sudo docker container run hello-world


sudo docker container ls

sudo docker container run --publish 80:80 -detach --name wechat nginx 
sudo docker container run --publish 8080:80 -detach --name mouse httpd 
sudo docker container stop mouse ; sudo docker container rm mouse 
sudo docker container stop wechat ; sudo docker container rm wechat
sudo docker container stop database ; sudo docker container rm database
sudo docker container run -d --publish 3306:3306 --name database -e MYSQL_RANDOM_ROOT_PASSWORD=yes mysql
sudo docker container logs database
# GENERATED ROOT PASSWORD: ainou9yooQu1Raefo4ijiseeQu4tahri

sudo docker container stop <image-id>
sudo docker container rm <image-id>
sudo docker container stats 
sudo docker container top <image-id>
sudo docker container inspect <image-id>

sudo docker container run -it 
sudo docker container start -ai wechat
sudo docker container exec -it wechat bash
sudo docker pull alpine

sudo docker container inspect --format '{{ .NetworkSettings.IPAddress }}' wechat

sudo docker network ls
sudo docker network inspect bridge
sudo docker network create my_app_net
sudo docker network ls
sudo docker network create --bridge my_app_net
sudo docker network connect
sudo docker network disconnect


# centos:7
# ubuntu:14.04
sudo docker container run -detach --name my_centos centos:7
sudo docker container run -detach --name my_ubuntu ubuntu:14.04


sudo docker container exec -it my_centos bash
yum update curl
curl --version
sudo docker container exec -it my_ubuntu bash
apt-get update && apt-get install curl
sudo docker container stop my_centos ; sudo docker container rm my_centos
sudo docker container stop my_ubuntu ; sudo docker container rm my_ubuntu