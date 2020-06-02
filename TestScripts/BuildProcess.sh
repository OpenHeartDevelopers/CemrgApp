#!/bin/sh

#Prepare
#sudo add-apt-repository --yes ppa:ubuntu-sdk-team/ppa
#sudo apt-get update -qq
#sudo apt-get build-dep qt5-default
#sudo apt-get install libxcb-xinerama0-dev

#Build QT
git clone --branch 5.12 https://code.qt.io/qt/qt5.git
cd qt5; git submodule update --init --recursive
../qt5/configure -developer-build -opensource -nomake examples -nomake tests
make -j$(nproc)

echo $PWD
ls

#Build MITK
git clone --branch v2018.04.2 https://phabricator.mitk.org/source/mitk.git


#mkdir Build; cd Build
#cmake ../MITK
#make -j$(nproc); cd ..
#cmake -DMITK_EXTENSION_DIRS:PATH=??? ../
