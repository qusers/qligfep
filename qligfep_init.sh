#!/bin/bash

export QLIGFEP="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export PATH=${QLIGFEP}:${PATH}

export PYTHONPATH=${QLIGFEP}:${PYTHONPATH}


echo "Please enter the absolute path to your Q directory (e.g., /home/user/software/q6):"
read Q_PATH

echo "Please enter the absolute path to your Schrödinger directory (e.g., /home/user/apps/schrodinger2024-4/):"
read SCHROD

echo "Please enter the name of your default HPC cluster (e.g., MYCLUSTER):"
read CLUSTER_NAME

SETTINGS_PY="${QLIGFEP}/settings.py"

sed -i "s|\\\$Q_PATH|${Q_PATH}|g" "$SETTINGS_PY"
sed -i "s/\$DEFAULT/${CLUSTER_NAME}/g" "$SETTINGS_PY"
sed -i "s|\\\$SCHROD|${SCHROD}|g" "$SETTINGS_PY"

echo
echo "settings.py updated with your Q directory, Schrödinger directory and default HPC cluster name."
echo "Please check setting.py to add more HPC clusters and cluster options"
