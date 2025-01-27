#!/bin/bash

export QLIGFEP="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export PATH=${QLIGFEP}:${PATH}

export PYTHONPATH=${QLIGFEP}:${PYTHONPATH}
