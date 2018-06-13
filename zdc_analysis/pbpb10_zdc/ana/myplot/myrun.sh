rm plots/*.pdf
rm nohup.out
nohup root -l -b plotAll.C+ > nohup.out &
