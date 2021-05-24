#!/usr/bin/env bash
OUTFILE="./CemrgApp/Applications/MainApp/version.txt"

MYVAR=$(git tag)
echo $MYVAR | awk -F v '{print $NF}' > $OUTFILE

MYVAR2=$(echo $MYVAR | awk -F v '{print $NF}')
git log v$MYVAR2 --oneline -1 | awk -F ' ' '{print $1}' >> $OUTFILE
