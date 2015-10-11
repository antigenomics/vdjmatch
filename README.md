# Alpha version of the documentation

# VDJdb: a curated databased of annotated V(D)J junctions

## Installing VDJdb

### Using Docker Image

[Docker](https://www.docker.com/) is an open platform for building, shipping and running distributed applications. It gives programmers, development teams and operations engineers the common toolbox they need to take advantage of the distributed and networked nature of modern applications.

- [Install Docker](https://docs.docker.com/installation/)
- Download and unrar any latest compiled docker release of VDJdb on [GitHub](https://github.com/antigenomics/vdjdb/releases)
- Compile docker image using next command: `docker build -t antigenomics/vdjdb .`
- Run docker image using `docker run `docker run -t -u antigenomics/vdjdb

### Compiling from sources

Compilation is performed via [Apache Maven](http://maven.apache.org/):

```bash
git clone https://github.com/mikessh/vdjdb.git
cd vdjdb
mvn clean install
```

Also you should install PostgreSQL server on your computer and run `dump.sql` script using your username and password.
If you don't know how to install and configure PostgreSQL server please see [corresponding PostgreSQL documentation section](https://wiki.postgresql.org/wiki/Detailed_installation_guides)

## Usage

### Docker image

General way to execute VDJdb using docker image would be the following,

```bash
docker run -t -i antigenomics/vdjdb
```

### Using command line

```bash
java -jar vdjdb-X.X-X.jar [OPTIONS]
```

## Options

| Shorthand | Long name       | Required                                     | Description                                                |
|-----------|-----------------|----------------------------------------------|------------------------------------------------------------|
| -h        |                 | No                                           | Display help message                                       |
| -d        | --database      | Yes (only if you are not using docker image) | PostgreSQL database name                                   |
| -u        | --user          | Yes (only if you are not using docker image) | PostgreSQL user name                                       |
| -p        | --password      | Yes (only if you are not using docker image) | PostgreSQL user password                                   |
| -e        | --errors        | No                                           | Show detailed information about errors                     |
| -sf       | --showFields    | No                                           | Display the available fields for which you can use filters |
| -fm       | --filterMatch   | No                                           | Match Filters                                               |
| -fp       | --filterPattern | No                                           | Pattern Filters                                             |
| -ff       | --fuzzyFilter   | No                                           | Fuzzy Filters                                               |


## Filters

### Match Filter

Basic usage: `-fm fieldName=string`\
For example: `-fm V=TRBV29-1`\
Also you can specify the third additional parameter: `-fm fieldName=string=false`. This filter will find all records in database which not match `-fm fieldName=string` filter.

### Pattern Filter

Basic usage: `-fp fieldName=pattern`\
For example: `-fp V=TRBV2_-1`\
Also you can specify the third additional parameter: `-fp fieldName=pattern=false`. This filter will find all records in database which not match `-fp fieldName=pattern` filter.

### Fuzzy Filter

Basic usage: `-ff fieldName=string=distance`\
For example: `-ff cdr3=CSARGQGRDAAF=3`




