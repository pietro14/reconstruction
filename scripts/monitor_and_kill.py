#!/bin/env python
from time import time, sleep
import subprocess,os,sys

if __name__ == "__main__":

    seconds = 900
    myPID = os.getpid()
    
    while True:
        sleep(seconds)
        cmd = "top -b -n 1 | grep python | grep $USER | grep -v {mypid}".format(mypid=myPID)
        result = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
        totCpu=0
        for line in result.splitlines():
            totCpu = totCpu + float(line.split()[8])
        if totCpu<0.1:
            print ("KILLING ALL PYTHON3.8 COMMANDS!")
            killCmd = 'killall -9 python3.8'
            os.system(killCmd)
            break
    sys.exit()

        
