#!/usr/bin/python
from os import system as sys

print "Hello"

Runs=[286665,286711,286717,286834,286854,286908,286967,286990,286995,287038,287044,287068,287222,287224,287259,287270,287281,287321,287330,287334,287378,287380,287382,287560,287594,287632,287706,287728,287827,287843,287866,287924,287931]
Runs=[286665,286711,286717,286834,286854,286908,286967,286990,286995,287038]

#MaxFiles=30
MaxFiles=1
EvPerRun=10000000
EvPerJob=200000
NJobPerRun=EvPerRun/EvPerJob
print "Jobs per file = %d" %NJobPerRun
NDiv=1

'''
#Test
MaxFiles=2
NDiv=1
'''

cnt=0
cnt_run=0
for run in Runs:
    cnt_run+=1
    for iev in range(0,NJobPerRun):
        sys("mkdir ./log/log_%d_%d" %(cnt_run,iev))
        Logfile="./log/log_%d_%d" %(cnt_run,iev)
        Filenm01="%s/AB_job_%d_%d.job" %(Logfile,cnt_run,iev)
        Filenm02="%s/AB_job_%d_%d.sh" %(Logfile,cnt_run,iev)
        
        sys("cp ~/condor_script4 %s" %Filenm01)
        sys("echo \"Output          = %s/Corrjob.out\">>%s" %(Logfile,Filenm01) )
        sys("echo \"Error           = %s/Corrjob.err\">>%s" %(Logfile,Filenm01) )
        sys("echo \"Log             = %s/Corrjob.log\">>%s" %(Logfile,Filenm01) )
        sys("echo \"Executable      = %s\">>%s" %(Filenm02,Filenm01) )
        sys("echo \"Queue\"                        >>%s" %Filenm01)
        
        sys("echo \"#!/bin/sh    \">%s" %Filenm02)
        sys("echo 'echo $0'>>%s" %Filenm02)
        sys("echo \"sh myrun_res.sh %d %d > %s/job.out\">>%s" %(run,iev,Logfile,Filenm02) )
        sys("echo \"rm %s\">>%s" %(Filenm01,Filenm02))
        sys("echo \"rm %s\">>%s" %(Filenm02,Filenm02))
        
        sys("chmod 711 %s" %Filenm02)
        sys("condor_submit %s" %Filenm01) 
        cnt+=1

print "All jobs submitted : %d" %cnt

'''
for fr in range (0,MaxFiles,NDiv):
    to=fr+NDiv
    for iev in range(0,Njobperfile):
        #if(cnt>5) : break
        #sys("rm ./output/output_%d_%d" %(fr,to))
        #sys("rm -rf ./log/log_%d_%d" %(fr,to))
        #print "We now delete log"
        sys("mkdir ./log/log_%d_%d_%d" %(fr,to,iev))
        Logfile="./log/log_%d_%d_%d" %(fr,to,iev)
        Filenm01="%s/FBasy_%d_%d_%d.job" %(Logfile,fr,to,iev)
        Filenm02="%s/FBasy_%d_%d_%d.sh" %(Logfile,fr,to,iev)
        
        sys("cp ~/condor_script4 %s" %Filenm01)
        sys("echo \"Output          = %s/Corrjob.out\">>%s" %(Logfile,Filenm01) )
        sys("echo \"Error           = %s/Corrjob.err\">>%s" %(Logfile,Filenm01) )
        sys("echo \"Log             = %s/Corrjob.log\">>%s" %(Logfile,Filenm01) )
        sys("echo \"Executable      = %s\">>%s" %(Filenm02,Filenm01) )
        sys("echo \"Queue\"                        >>%s" %Filenm01)
        
        sys("echo \"#!/bin/sh    \">%s" %Filenm02)
        sys("echo 'echo $0'>>%s" %Filenm02)
        sys("echo \"sh myrun_res.sh %d %d %d > %s/job.out\">>%s" %(fr,to,iev,Logfile,Filenm02) )
        sys("echo \"rm %s\">>%s" %(Filenm01,Filenm02))
        sys("echo \"rm %s\">>%s" %(Filenm02,Filenm02))
        
        sys("chmod 711 %s" %Filenm02)
        sys("condor_submit %s" %Filenm01) 
        cnt+=1

print "All jobs submitted : %d" %cnt
'''
