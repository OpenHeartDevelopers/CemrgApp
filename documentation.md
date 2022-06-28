---
layout: page
title: Documentation
permalink: /documentation/
---
**Contents:**
- [Installation](#installation)
    - [Troubleshooting linux](#troubleshooting-linux)
- [Basic Usage](#basic-usage)
- [Advanced CemrgApp functionalities through Docker](#advanced-cemrgapp-functionalities-through-docker)

## Installation 
CemrgApp comes packaged in binaries, which can be installed easily.
Follow the instruction depending on your operating system:
+ **Windows:** Unzip the file to a folder (e.g `C:\Users\YOUR_NAME\apps`) and open the `CemrgApp.exe` file. You can right-click and select Add to taskbar to make it easily accessible. 
+ **macOS:** Open the DMG file, then drag CemrgApp to your Applications folder. 
+ **linux:** Unzip the file to a folder (e.g `~/apps`) and run the `CemrgApp.sh` on a terminal.

#### Troubleshooting linux
You might need to install some libraries to run CemrgApp appropriately, 
if you run into some problem, consider downloading some of the following libraries
```sh
sudo apt install lib-tbb
```
> Notice: these instructions are only useful in debian-based distributions, check your distribution if you're unsure

## Basic Usage 
CemrgApp has different perspectives, small apps that do one functionality. 
Check the Examples page to learn more. Currently, the available perspectives are: 

+ Motion Quantification
+ Anatomical Measurements
+ Scar Quantification
+ Morphological Measurements
+ Electrophysiology Simulations

## Advanced CemrgApp functionalities through Docker 

Installing Docker is necessary to use CemrgApp's Advanced Modules, 
like the automatic segmentation in the Scar Plugin. 
Docker can be downloaded or installed from:

- [Docker](https://docs.docker.com/install/linux/docker-ce/ubuntu) for Linux.
  - You need to provide access to Docker by following these [instructions](https://docs.docker.com/install/linux/linux-postinstall), or see below for a summarised version of the instructions.
  - WARNING: this will give root access to Docker containers.
- [Docker](https://docs.docker.com/docker-for-mac/install) for macOS.
  - Simple installation and running.
- [Docker](https://docs.docker.com/docker-for-windows/install) for Microsoft Windows.

After installing, you will need to configure docker. 
The best guide is in our [wiki page](https://github.com/CemrgAppDevelopers/CemrgApp/wiki/Configure-Docker-with-CemrgApp). 