inputpath="../inputlist/rootfiles_1.list"
root -l <<EOF
   .L extractor_flat.C+
   extractor_flat* extr1 = new extractor_flat("../inputlist/rootfiles_1.list");
   extr1.exec();
   .q
EOF


