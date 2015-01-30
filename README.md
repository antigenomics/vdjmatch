# VDJdb: a curated databased of annotated V(D)J junctions

## Purpose

This repository contains the database table itself stored as a plain-text file. This is done in order to manage database versions via git.
The repository also contains API for database browsing, see [Vdjtools](https://github.com/mikessh/vdjtools) for performing batch queries for immune repertoire samples.

## Installation

The only required dependency is [MiLib](https://github.com/milaboratory/milib).
Installation is performed via [Apache Maven](http://maven.apache.org/). Should be installed under Java v1.7.

```bash
git clone --branch 1.0 --depth 1 https://github.com/milaboratory/milib.git
cd milib && mvn clean install && cd ..
git clone https://github.com/mikessh/vdjdb.git
cd vdjdb && mvn clean install
```

