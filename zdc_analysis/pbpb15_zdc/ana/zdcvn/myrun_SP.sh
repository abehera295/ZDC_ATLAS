#!/bin/bash
I1=$1
I2=$2
IEv=$3
NEv=200000

root -b <<EOF
.L extractor_flat.C+g
gSystem->Load("extractor_flat_C.so")
extractor_flat *extr1 = new extractor_flat("../inputlist/rootfiles.list",$I1,$I2,$IEv,$NEv);
extr1->Exec_SP();
.q
EOF

#gSystem->Load("extractor_flat_C.so");
