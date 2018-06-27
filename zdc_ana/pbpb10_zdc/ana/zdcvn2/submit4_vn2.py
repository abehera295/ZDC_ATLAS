#!/usr/bin/python
from os import system as sys

print "Hello"


MaxFiles=30
#MaxFiles=1
Evtsperfile=3200000
Evtsperjob=100000
Njobperfile=Evtsperfile/Evtsperjob
print "Jobs per file = %d" %Njobperfile
NDiv=1

'''
#Test
MaxFiles=2
NDiv=1
'''

cnt=0
for fr in range (0,MaxFiles,NDiv):
    to=fr+NDiv
    for iev in range(0,Njobperfile):
        #sys("rm ./output/output_%d_%d" %(fr,to))
        #sys("rm -rf ./log/log_%d_%d" %(fr,to))
        #print "We now delete log"
        sys("mkdir ./log/log_%d_%d_%d" %(fr,to,iev))
        Logfile="./log/log_%d_%d_%d" %(fr,to,iev)
        Filenm01="%s/FBasy_%d_%d_%d.job" %(Logfile,fr,to,iev)
        Filenm02="%s/FBasy_%d_%d_%d.sh" %(Logfile,fr,to,iev)
    
        sys("cp ~/condor_script2 %s" %Filenm01)
        sys("echo \"Output          = %s/Corrjob.out\">>%s" %(Logfile,Filenm01) )
        sys("echo \"Error           = %s/Corrjob.err\">>%s" %(Logfile,Filenm01) )
        sys("echo \"Log             = %s/Corrjob.log\">>%s" %(Logfile,Filenm01) )
        sys("echo \"Executable      = %s\">>%s" %(Filenm02,Filenm01) )
        sys("echo \"queue\"                        >>%s" %Filenm01)
        
        sys("echo \"#!/bin/sh    \">%s" %Filenm02)
        sys("echo 'echo $0'>>%s" %Filenm02)
        sys("echo \"sh myrun_vn2.sh %d %d %d > %s/job.out\">>%s" %(fr,to,iev,Logfile,Filenm02) )
        sys("echo \"rm %s\">>%s" %(Filenm01,Filenm02))
        sys("echo \"rm %s\">>%s" %(Filenm02,Filenm02))
        
        sys("chmod 711 %s" %Filenm02)
        sys("condor_submit %s" %Filenm01) 
        cnt+=1

print "All jobs submitted : %d" %cnt

