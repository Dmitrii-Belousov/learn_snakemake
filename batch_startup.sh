#!/bin/bash
# Azure Batch startup script to install Apptainer

set -e

echo "Installing Apptainer..."

# Update package lists
apt-get update

# Install dependencies
apt-get install -y \
    wget \
    build-essential \
    libseccomp-dev \
    pkg-config \
    uidmap \
    squashfs-tools \
    fakeroot \
    cryptsetup \
    dh-apparmor

# Install Go
cd /tmp
wget https://go.dev/dl/go1.21.5.linux-amd64.tar.gz
rm -rf /usr/local/go
tar -C /usr/local -xzf go1.21.5.linux-amd64.tar.gz
export PATH=$PATH:/usr/local/go/bin

# Install Apptainer
cd /tmp
wget https://github.com/apptainer/apptainer/releases/download/v1.3.0/apptainer-1.3.0.tar.gz
tar -xzf apptainer-1.3.0.tar.gz
cd apptainer-1.3.0

./mconfig
make -C builddir
make -C builddir install

echo "Apptainer installation complete"
apptainer --version
