#!/bin/sh

#Build QT
sudo apt-get build-dep qt5-default
sudo apt-get install libxcb-xinerama0-dev

#Build MITK
git clone --branch v2018.04.2 https://phabricator.mitk.org/source/mitk.git MITK

#Download CemrgApp
#git


mkdir Build; cd Build
#cmake -DMITK_EXTENSION_DIRS:PATH=??? ../
cmake ../MITK
make -j8
