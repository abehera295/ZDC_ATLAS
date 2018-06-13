root -l <<EOF
   .L extractor_flat.C+
   extractor_flat* extr1 = new extractor_flat("input.txt");
   extr1.exec();
   .q
EOF


