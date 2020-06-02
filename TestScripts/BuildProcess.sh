#!/bin/sh


#Build MITK
git clone --branch v2018.04.2 https://phabricator.mitk.org/source/mitk.git MITK
mkdir Build; cd Build
cmake ../MITK
make -j$(nproc)
#cmake -DMITK_EXTENSION_DIRS:PATH=??? ../
