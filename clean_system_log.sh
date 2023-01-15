#! /bin/bash
 
# clean journals
rm -rf /var/log/journal/*
 
# clean .gz logs
# ls -hlS /var/log/ | grep -E '.gz'
ls /var/log/ | grep -E '.gz' | xargs -i rm -r /var/log/{}
 
# Echo nothing to logs
# 这里并未将所有的log去除，只是将那些体积比较大的log去除。强迫症患者可以自己添加一些项目进去，比如tallylog/faillog之类的
ls /var/log/ | grep -E 'syslog|messages|user|daemon|btmp|auth|mail' | xargs -i tee /var/log/{}
 
# Report
echo 'Clean all system logs!'
