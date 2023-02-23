#!/usr/bin/env bash
OUTFILE="./CemrgApp/Applications/MainApp/version.txt"

# SHAVAR=$(git log --oneline -1 | awk -F ' ' '{print $1}')
# VERVAR=$(git tag --points-at=${SHAVAR})
#
# echo $VERVAR
# echo $SHAVAR
#
# echo $VERVAR > $OUTFILE
# echo $SHAVAR >> $OUTFILE

MYVAR=$(git tag)
echo $MYVAR | awk -F v '{print $NF}' > $OUTFILE

MYVAR2=$(echo $MYVAR | awk -F v '{print $NF}')
git log v$MYVAR2 --oneline -1 | awk -F ' ' '{print $1}' >> $OUTFILE
