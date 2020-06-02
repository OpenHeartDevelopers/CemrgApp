#!/bin/sh

#Build QT
apt-get build-dep qt5-default
apt-get install libxcb-xinerama0-dev

git clone https://phabricator.mitk.org/source/mitk.git MITK
git checkout v2018.04.2

#Download CemrgApp
#git


mkdir Build; cd Build
#cmake -DMITK_EXTENSION_DIRS:PATH=??? ../
cmake ..
make -j8
