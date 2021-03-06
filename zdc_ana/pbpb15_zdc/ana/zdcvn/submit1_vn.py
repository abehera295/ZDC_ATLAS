#!/usr/bin/python
from os import system as sys

print "Hello"

Runs=[286665,286711,286717,286834,286854,286908,286967,286990,286995,287038,287044,287068,287222,287224,287259,287270,287281,287321,287330,287334,287378,287380,287382,287560,287594,287632,287706,287728,287827,287843,287866,287924,287931]
Runs=[286665,286711,286717,286834,286854,286908,286967,286990,286995,287038]
NEvts=[200000,1400000,3600000,3800000,8400000,5200000,5800000,4800000,800000,6200000]
#Runs=[286711]

#MaxFiles=30
#MaxFiles=1
#EvPerRun=2000000
EvPerJob=100000
#NJobPerRun=EvPerRun/EvPerJob
#print "Jobs per file = %d" %NJobPerRun
NDiv=1

'''
#Test
MaxFiles=2
NDiv=1
'''

cnt=0
cnt_run=0
for run in Runs:
    NEv=NEvts[Runs.index(run)]
    NJobs=NEv/EvPerJob
    cnt_run+=1
    for iev in range(0,NJobs):
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
        sys("echo \"sh myrun_vn.sh %d %d %d > %s/job.out\">>%s" %(run,iev,EvPerJob,Logfile,Filenm02) )
        sys("echo \"rm %s\">>%s" %(Filenm01,Filenm02))
        sys("echo \"rm %s\">>%s" %(Filenm02,Filenm02))
        
        sys("chmod 711 %s" %Filenm02)
        sys("condor_submit %s" %Filenm01) 
        cnt+=1

print "All jobs submitted : %d" %cnt

