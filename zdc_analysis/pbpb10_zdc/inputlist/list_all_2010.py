#!/usr/bin/python
import os
Run=1
dirin=[]
dirin.append("/Users/arabinda/StonyBrook/zdc_v1/output_DS/user.arbehera.pDST.28Aug2016.169927_PDST.root")
dirin.append("/Users/arabinda/StonyBrook/zdc_v1/output_DS/user.arbehera.pDST.21Dec2017.169927_PDST.root")
os.system("rm rootfiles*")

for runpath in dirin:
    runno=dirin.index(runpath)+1
    cmd=("ls %s/*.root > rootfiles_%d.list"%(runpath,runno))
    os.system("%s"%cmd)


