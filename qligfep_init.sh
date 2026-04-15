#!/bin/bash

export QLIGFEP="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export PATH=${QLIGFEP}:${PATH}

export PYTHONPATH=${QLIGFEP}:${PYTHONPATH}

echo "========================================="
echo "QligFEP Installation Setup"
echo "========================================="
echo

# Check if Q is available locally in the repository
LOCAL_Q_PATH="${QLIGFEP}/Q/bin"
if [ -f "${LOCAL_Q_PATH}/qprep" ] || [ -f "${LOCAL_Q_PATH}/qdyn" ]; then
    echo "Found Q installation in local repository: ${LOCAL_Q_PATH}"
    read -p "Use this Q installation? (y/n) [y]: " USE_LOCAL
    USE_LOCAL=${USE_LOCAL:-y}
    if [[ "$USE_LOCAL" == "y" ]]; then
        Q_PATH="${LOCAL_Q_PATH}"
    else
        read -p "Enter the absolute path to your Q directory (e.g., /home/user/software/q6/bin): " Q_PATH
    fi
else
    read -p "Enter the absolute path to your Q directory (e.g., /home/user/software/q6/bin): " Q_PATH
fi

# Validate Q installation
if [ ! -f "${Q_PATH}/qprep" ] && [ ! -f "${Q_PATH}/qdyn" ]; then
    echo "ERROR: Q binaries not found in ${Q_PATH}"
    echo "Please ensure qprep and/or qdyn are available in the specified directory."
    exit 1
fi

read -p "Enter the absolute path to your Schrödinger directory (leave blank if not installed): " SCHROD_DIR
if [ -z "${SCHROD_DIR}" ]; then
    SCHROD_DIR=""
fi

read -p "Enter the name of your default HPC cluster (e.g., MYCLUSTER) [default]: " CLUSTER_NAME
CLUSTER_NAME=${CLUSTER_NAME:-default}

# Replace placeholders in settings.py
SETTINGS_PY="${QLIGFEP}/settings.py"
TEMPLATE_FILE="${QLIGFEP}/settings.py.template"

# Restore template if it exists
if [ -f "${TEMPLATE_FILE}" ]; then
    cp "${TEMPLATE_FILE}" "${SETTINGS_PY}"
else
    echo "WARNING: Template file not found."
fi

# Use different delimiters for sed to handle paths with slashes
sed -i "s|\$Q_PATH|${Q_PATH}|g" "$SETTINGS_PY"
sed -i "s|\$SCHROD_DIR|${SCHROD_DIR}|g" "$SETTINGS_PY"
sed -i "s|\$CLUSTER_NAME|${CLUSTER_NAME}|g" "$SETTINGS_PY"

echo
echo "========================================="
echo "Setup Complete!"
echo "========================================="
echo "Q directory: ${Q_PATH}"
echo "Schrödinger directory: ${SCHROD_DIR:-not installed}"
echo "Default HPC cluster: ${CLUSTER_NAME}"
echo
echo "You can now run QligFEP, QresFEP, or QLIE commands."
echo "To add more HPC clusters, edit settings.py manually."
