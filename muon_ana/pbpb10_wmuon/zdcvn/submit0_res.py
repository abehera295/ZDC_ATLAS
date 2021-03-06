#!/usr/bin/python
from os import system as sys

print "Hello"

#2010 runs
Runs=[170482,170459,170398,170082,170080,170016,170015,170004,170002,169966,169964,169961,169927,169864,169839,169783,169765,169751,169750,169693,169627,169567,169566,169564,169270,169226,169224,169223,169207,169206,169175,169136]
#2015 runs
#Runs=[286665,286711,286717,286834,286854,286908,286967,286990,286995,287038,287044,287068,287222,287224,287259,287270,287281,287321,287330,287334,287378,287380,287382,287560,287594,287632,287706,287728,287827,287843,287866,287924,287931]
skipRuns=[]

NEvts=[3000000]
#Runs=[286711]

#MaxFiles=30
#MaxFiles=1
EvPerRun=3000000
EvPerJob=200000
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
    if(run in skipRuns): continue
    #NEv=NEvts[Runs.index(run)]
    NJobs=EvPerRun/EvPerJob
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
        sys("echo \"sh myrun_res.sh %d %d %d > %s/job.out\">>%s" %(run,iev,EvPerJob,Logfile,Filenm02) )
        sys("echo \"rm %s\">>%s" %(Filenm01,Filenm02))
        sys("echo \"rm %s\">>%s" %(Filenm02,Filenm02))
        
        sys("chmod 711 %s" %Filenm02)
        sys("condor_submit %s" %Filenm01) 
        cnt+=1

print "All jobs submitted : %d" %cnt

