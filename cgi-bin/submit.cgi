#!/usr/bin/python -u

import os
import sys
import time
import uuid
import shlex
import subprocess
import traceback

import cgi
import cgitb; cgitb.enable()

#
if __name__ == '__main__':
    print "Content-Type: text/html"
    print
    
    try:
        HTMLForm    = cgi.FieldStorage()
        HTMLTime    = time.gmtime()
        HTMLSession = None
        
        while True:
            HTMLSession = str(uuid.uuid4())
            if not os.path.exists("../tmp/" + HTMLSession + ".txt"):
                break
        
        Data = []
        Data.append("# IP: " + str(os.environ['REMOTE_ADDR']))
        Data.append("# Date: " + time.strftime("%d-%m-%Y", HTMLTime))
        Data.append("# Time: " + time.strftime("%H:%M:%S", HTMLTime))
        Data.append("# Weights: " + str(HTMLForm['weighted'].value))
        
        if HTMLForm.has_key("threshold") and HTMLForm['threshold'].value:
            Data.append("# Threshold: " + str(HTMLForm['threshold'].value))
        if HTMLForm.has_key("phenotypes") and HTMLForm['phenotypes'].value:
            Data.append("# Phenotypes: " + str(HTMLForm['phenotypes'].value))
        
        Data.append(str(HTMLForm['batch'].value.strip()))
        
        open("../tmp/" + HTMLSession + ".txt", "w").write(
            "\n".join(Data)
        )
        
        run      = "python fathmm.py ../tmp/{0}.txt /tmp/{0}.tab -w {1} ".format(HTMLSession, str(HTMLForm['weighted'].value))
        if HTMLForm.has_key("threshold") and HTMLForm['threshold'].value:
            run += "-t {0} ".format(str(HTMLForm['threshold'].value))
        if HTMLForm.has_key("phenotypes") and HTMLForm['phenotypes'].value:
            run += "-p {0} ".format(str(HTMLForm['phenotypes'].value))
            
        mv       = "mv /tmp/{0}.tab ../tmp/{0}.tab".format(HTMLSession)
        
        #print str(HTMLForm) + "<br \>"
        #print run + "<br \>"
        #print mv
        
        pid = subprocess.Popen("{0}; {1}".format(run, mv),
                               shell     = True,
                               stdin     = subprocess.PIPE, 
                               stdout    = subprocess.PIPE,
                               stderr    = subprocess.STDOUT,
                               close_fds = True
                              )
        
        print """
<html>
    <head>
        <meta charset="utf-8">
        <title>fathmm - New Submission</title>
        <script type="text/JavaScript">
            function redirect() {
                window.location = './results.cgi?session=""" + HTMLSession + """';
            }
        </script>
    </head>
    <body onLoad="setTimeout('redirect()', 0)">
        &nbsp;
    </body>
</html>
    """
    except Exception, e:
        raise
        
