#!/bin/bash

# Check if Go is already installed
check_go_installed() {
  if command -v go &> /dev/null; then
    echo "Go is already installed."
    go version
    exit 0
  else
    echo "Go is not installed. Installing the latest version..."
  fi
}

# Get the latest version of Go
get_latest_go_version() {
  curl -s https://go.dev/VERSION?m=text
}

# Install Go on Linux
install_go_linux() {
  LATEST_VERSION=$(get_latest_go_version)
  echo "Latest Go version: $LATEST_VERSION"
  wget https://dl.google.com/go/${LATEST_VERSION}.linux-amd64.tar.gz
  sudo tar -C /usr/local -xzf ${LATEST_VERSION}.linux-amd64.tar.gz
  rm ${LATEST_VERSION}.linux-amd64.tar.gz
  echo "export PATH=\$PATH:/usr/local/go/bin" >> ~/.profile
  source ~/.profile
}

# Install Go on macOS
install_go_macos() {
  LATEST_VERSION=$(get_latest_go_version)
  echo "Latest Go version: $LATEST_VERSION"
  wget https://dl.google.com/go/${LATEST_VERSION}.darwin-amd64.tar.gz
  sudo tar -C /usr/local -xzf ${LATEST_VERSION}.darwin-amd64.tar.gz
  rm ${LATEST_VERSION}.darwin-amd64.tar.gz
  echo "export PATH=\$PATH:/usr/local/go/bin" >> ~/.bash_profile
  source ~/.bash_profile
}

# Detect the operating system
install_go() {
  OS=$(uname -s)

  case "$OS" in
    Linux*)
      install_go_linux
      ;;
    Darwin*)
      install_go_macos
      ;;
    *)
      echo "Unsupported OS: $OS"
      exit 1
      ;;
  esac

  # Check if Go was installed correctly
  if command -v go &> /dev/null; then
    echo "Go successfully installed!"
    go version
  else
    echo "Failed to install Go."
    exit 1
  fi
}

# Main execution starts here
check_go_installed
install_go
go get -u ./...