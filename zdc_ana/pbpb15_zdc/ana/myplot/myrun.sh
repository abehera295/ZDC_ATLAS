rm plots/*.pdf
rm nohup.out
nohup root -l -b plotall.C+ > nohup.out &
