I1=0
I2=0
IEv=0
NEv=200000
root -b <<EOF
.L extractor_flat.C+
extractor_flat *extr1 = new extractor_flat("../inputlist/rootfiles.list",$I1,$I2,$IEv,$NEv);
extr1->Exec_doAna();
.q
EOF

#gSystem->Load("extractor_flat_C.so");
