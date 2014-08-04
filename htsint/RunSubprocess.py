#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,subprocess, threading

def run_subprocess(cmd):
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    while True:
        try:
            next_line = proc.stdout.readline()
            proc.wait()
            if next_line == '' and proc.poll() != None:
                break
        except:
            proc.wait()
            break


class RunSubprocess(object):
    """
    a generic class 
    """

    def __init__(self, cmd, mainWindow=None):
        self.cmd = cmd
        self.mainWindow = mainWindow
        self.process = None
        self.stdOut,self.stdErr = None,None

    def run(self,timeout=100):
        def target():
            self.process = subprocess.Popen(self.cmd,shell=True,stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE,universal_newlines=True,bufsize=4096)

            self.stdOut, self.stdErr = self.process.communicate()

        self.thread = threading.Thread(target=target)
        self.thread.start()

        ## wait a specified amount of time before terminating
        if timeout != None:
            self.thread.join(timeout)
            if self.thread.is_alive():
                print 'The subprocess was auto-terminated due to timeout'
                print "...", self.process.poll()
                self.process.terminate()
                self.thread.join()
        
            return self.process.returncode
        return None

    def add_text(self):
        if self.textScreen == None:
            return
        
        if self.textScreen.showMessages == False:
            self.textScreen.toggle_message_btn()
            
        self.textScreen.add_text(txt)
        QtGui.QApplication.processEvents()

    def terminate(self):
        if self.thread.is_alive():
            self.process.terminate()
            self.thread.join()

if __name__ == '__main__':
    
    myProcess = RunSubprocess("echo 'Process started'; sleep 2; echo 'Process finished'")
    
    ## test should pass
    returnCode = myProcess.run(timeout=10)
    print 'pass return code', returnCode

    ## test should fail
    returnCode = myProcess.run(timeout=1)
    print 'fail return code', returnCode
