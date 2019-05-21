# Analysis of spectral search hits

This directory contains the script used to analyze the output of Ionbot search tool

There is a jupyter notebook to test out functions and figures on a subset of the output. Once successful they are moved into ionbot_analysis.py for use on all the data.

I run the jupyter notebook on a docker, which I bring up with the following commands:

docker-compose build
docker-compose up -d
docker exec -it scripts_jupyter_1 bash -c 'jupyter notebook list' | grep http | cut -f1 -d ' ' | awk '{sub(/8888/, "8887")}1' | xargs open
