#!/usr/bin/python

import os

print "Enter choice :"
print "Choice = 1 --> ZDC Flattening,"
print "Choice = 2 --> Eta-Phi Wights,"
print "Choice = 3 --> v1 Calculation"
n=input()
myfile="./nohup.out"
if(os.path.isfile(myfile)==True):
    os.remove(myfile)
'''
if(int(n)==1):
    os.system("nohup sh myrun_zdc_psi.sh > nohup.out &")

if(int(n)==2):
    os.system("nohup sh myrun_wei.sh > nohup.out &")

if(int(n)==3):
    os.system("nohup sh myrun_v1.sh > nohup.out &")
'''

