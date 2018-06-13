RunNo=$1
IEv=$2
NEv=200000

root -b <<EOF
.L extractor_flat.C+g
gSystem->Load("extractor_flat_C.so")
extractor_flat *extr1 = new extractor_flat("../inputlist/rootfiles_$RunNo.list",$RunNo,$IEv,$NEv);
extr1->Exec_res();
.q
EOF

#gSystem->Load("extractor_flat_C.so");
