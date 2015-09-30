# vdjdb docker container. Version: 0.0.1
FROM ubuntu:14.04
MAINTAINER Dmitri Bagaev<bvdmitri@gmail.com>

# Install `add-apt-repository` and Java8 in a single command
RUN sudo apt-get update && \
	sudo apt-get install -y software-properties-common && \
	echo debconf shared/accepted-oracle-license-v1-1 select true | sudo debconf-set-selections && \
	echo debconf shared/accepted-oracle-license-v1-1 seen true | sudo debconf-set-selections && \
	sudo add-apt-repository -y ppa:webupd8team/java && \
	sudo apt-get update && \
	sudo apt-get install -y oracle-java8-installer && \
	sudo rm -rf /var/lib/apt/lists/* && \
	sudo rm -rf /var/cache/oracle-jdk8-installer


RUN mkdir /vdjdb
COPY vdjdb-1.0-SNAPSHOT.jar /vdjdb/
COPY dump.sql /vdjdb/
WORKDIR /vdjdb/

# Add the PostgreSQL PGP key to verify their Debian packages.
# It should be the same key as https://www.postgresql.org/media/keys/ACCC4CF8.asc
RUN apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys B97B0AFCAA1A47F044F244A07FCC7D46ACCC4CF8

# Add PostgreSQL's repository. It contains the most recent stable release
#     of PostgreSQL, ``9.3``.
RUN echo "deb http://apt.postgresql.org/pub/repos/apt/ precise-pgdg main" > /etc/apt/sources.list.d/pgdg.list

# Install ``python-software-properties``, ``software-properties-common`` and PostgreSQL 9.3
#  There are some warnings (in red) that show up during the build. You can hide
#  them by prefixing each apt-get statement with DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y python-software-properties software-properties-common postgresql-9.3 postgresql-client-9.3 postgresql-contrib-9.3

USER postgres

RUN    /etc/init.d/postgresql start &&\
    psql --command "CREATE USER docker WITH SUPERUSER PASSWORD 'docker';" &&\
    createdb -O docker docker &&\
    psql --dbname=docker --file=dump.sql

# Adjust PostgreSQL configuration so that remote connections to the
# database are possible. 
RUN echo "host all  all    0.0.0.0/0  md5" >> /etc/postgresql/9.3/main/pg_hba.conf

# And add ``listen_addresses`` to ``/etc/postgresql/9.3/main/postgresql.conf``
RUN echo "listen_addresses='*'" >> /etc/postgresql/9.3/main/postgresql.conf

# Expose the PostgreSQL port
EXPOSE 5432

# Add VOLUMEs to allow backup of config, logs and databases
VOLUME  ["/etc/postgresql", "/var/log/postgresql", "/var/lib/postgresql"]

COPY search.sh /vdjdb/
ENTRYPOINT ["sh", "search.sh"]
