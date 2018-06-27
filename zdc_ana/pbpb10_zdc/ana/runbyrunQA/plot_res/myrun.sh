rm nohup.out
rm plots_1/*.pdf
rm plots_2/*.pdf
nohup root -l -b draw_res.C+ > nohup.out &
