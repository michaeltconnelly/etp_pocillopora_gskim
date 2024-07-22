#!/bin/bash
# script for extracting the likelihood values for different Ks from NGSadmix log files

for log in `ls */*pgra*log`
do 
echo `basename $log` | sed 's/.*\(K[0-9]*\)_[0-9].*/\1/g'
grep "best" $log | sed 's/best like=\(-[0-9]*.[0-9]*\) after .* iterations/\1/g'
done | sed 'N;s/\n/\t/' | sed 's/^K//g' | sort -k 1,1n > ${set}_ngsadmix_likeKlogs.txt

