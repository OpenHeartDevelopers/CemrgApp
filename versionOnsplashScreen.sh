#!/usr/bin/env bash

MYVAR=$(git tag)
echo $MYVAR | awk -F v '{print $NF}' > ./CemrgApp/Applications/MainApp/version.txt
echo $MYVAR | awk -F v '{print $NF}'
