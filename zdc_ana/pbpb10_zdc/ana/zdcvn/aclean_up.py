#!/usr/bin/python
import os
from os import system as sys
import shutil

print "Do you want to delete rootfiles (yes/no)?"
choice=raw_input()
if (choice=="yes"):
    print "You choose yes. So rootfiles will be deleted."
elif (choice=="no"):
    print "You choose no. So rootfiles will not be deleted."
else:
    print "Wrong choice!"

if(os.path.exists("./nohup.out")): os.remove("./nohup.out")
os.system("rm tempout*")

if(os.path.exists("./log")): shutil.rmtree("./log")
if(os.path.exists("./log2")): shutil.rmtree("./log2")
if(os.path.exists("./rmlist")): shutil.rmtree("./rmlist")
os.makedirs("./log")
os.makedirs("./log2")
os.makedirs("./rmlist")

dir1="./out"

if(choice=='yes'):
    if(os.path.exists(dir1)): shutil.rmtree(dir1)
    os.makedirs(dir1)
    #os.makedirs(dir1+"/histoutput")
    #os.makedirs(dir1+"/treeoutput")
    

