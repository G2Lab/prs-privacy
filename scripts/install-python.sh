#!/bin/bash

install_python_linux() {
    echo "Detected Linux. Installing Python..."
    sudo apt update -y
    sudo apt install -y software-properties-common
    sudo add-apt-repository -y ppa:deadsnakes/ppa
    latest_python=$(apt-cache search "^python3\.[0-9]$" | awk '{print $1}' | sort -r | head -n1)
    sudo apt install -y "$latest_python"
    sudo ln -sf /usr/bin/"$latest_python" /usr/local/bin/python3
}

install_python_macos() {
    echo "Detected macOS. Installing Python..."
    brew install python
}

OS=$(uname)
if [[ "$OS" == "Linux" ]]; then
    install_python_linux
elif [[ "$OS" == "Darwin" ]]; then
    install_python_macos
else
    echo "Unsupported OS: $OS. Exiting..."
    exit 1
fi

if ! command -v pip3 &>/dev/null; then
    echo "Installing pip..."
    curl -sS https://bootstrap.pypa.io/get-pip.py | sudo -H python3
fi

echo "Creating virtual environment and installing required packages..."
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt