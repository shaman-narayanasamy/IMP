#!/bin/bash
set -e

IMP_VERSION=$1
docker build -t docker-r3lab.uni.lu/imp/imp-deps:$IMP_VERSION -f Dockerfile-dependencies .
docker build -t docker-r3lab.uni.lu/imp/imp-tools:$IMP_VERSION -f Dockerfile-tools .
docker build -t docker-r3lab.uni.lu/imp/imp:$IMP_VERSION .
docker save docker-r3lab.uni.lu/imp/imp:$IMP_VERSION > imp-$IMP_VERSION.tar
gzip imp-$IMP_VERSION.tar
