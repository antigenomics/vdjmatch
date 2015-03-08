# VDJdb: a curated databased of annotated V(D)J junctions

## Purpose

This repository contains the database table itself stored as a plain-text file. This is done in order to manage database versions via git.
The repository also contains API for database browsing, see [Vdjtools](https://github.com/mikessh/vdjtools) for performing batch queries for immune repertoire samples.
VDJdb extensively uses [MiLib](https://github.com/milaboratory/milib) NGS processing Java library. Should be installed and run under Java v1.7+.

## Compiling from sources

Compilation is performed via [Apache Maven](http://maven.apache.org/):

```bash
git clone https://github.com/mikessh/vdjdb.git
cd vdjdb
mvn clean install
```

