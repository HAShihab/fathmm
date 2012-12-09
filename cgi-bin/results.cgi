#!/usr/bin/python -u

import os
import sys
import time
import uuid
import shlex
import subprocess

import cgi
import cgitb; cgitb.enable()

#
if __name__ == '__main__':
    try:
        print "Content-Type: text/html"
        print
        
        HTMLForm    = cgi.FieldStorage()
        HTMLTime    = time.gmtime()
        HTMLSession = None
        
        if not HTMLForm.has_key('session'):
            while True:
                # assert a unique session is generated
                HTMLSession = str(uuid.uuid4())
                
                if not os.path.exists("../tmp/" + HTMLSession + ".txt"):
                    break
            
            # write the submission (& server parameters) to the server
            for x in [ "# Remote Address: " + str(os.environ['REMOTE_ADDR']) + "\n",
                       "# Server Date: " + time.strftime("%d-%m-%Y", HTMLTime) + "\n",
                       "# Server Date: " + time.strftime("%H:%M:%S", HTMLTime) + "\n",
                       str(HTMLForm['batch'].value)
                     ]:
                open("../tmp/" + HTMLSession + ".txt", "a").write(x)
            
            # perform our analysis/predictions in the background ...
            pid = subprocess.Popen([ sys.executable, 'fathmm.py', '-i', '../tmp/' + HTMLSession + '.txt', '-w', str(HTMLForm['weighted'].value), '-p', str(HTMLForm['phenotypes'].value), '-H' ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        else:
            HTMLSession = str(HTMLForm['session'].value)
        
        print """
<!DOCTYPE html>
<!--[if lt IE 7]> <html class="lt-ie9 lt-ie8 lt-ie7"> <![endif]-->
<!--[if IE 7]>    <html class="lt-ie9 lt-ie8"> <![endif]-->
<!--[if IE 8]>    <html class="lt-ie9"> <![endif]-->
<!--[if gt IE 8]><!--> <html class=""> <!--<![endif]-->
    <head>
        <title>fathmm - Inherited Mutation Analysis</title>

        <meta charset="utf-8">
        
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <meta name="description" content="Predict the Functional and Phenotypic Consequences of Protein Missense Variants">
        <meta name="keywords" content="Functional Analysis, Protein Missense Variants, Amino Acid Substitutions, Hidden Markov Models, HMMs, Single Nucleotide Polymorphisms, SNPs, Non Synonymous Mutation, nsSNPs, FATHMM">
        <meta name="author" content="Hashem A. Shihab">
        
        <meta http-equiv="cache-control" content="max-age=0" />
        <meta http-equiv="cache-control" content="no-cache" />
        <meta http-equiv="expires" content="0" />
        <meta http-equiv="expires" content="Tue, 01 Jan 1980 1:00:00 GMT" />
        <meta http-equiv="pragma" content="no-cache" />
        
        <!-- HTML5 "Upscaling" -->
        <!--[if lt IE 9]>
            <script src="http://html5shiv.googlecode.com/svn/trunk/html5.js"></script>
        <![endif]-->

        <link href="../css/bootstrap.css" rel="stylesheet">
        <link href="../css/bootstrap-responsive.css" rel="stylesheet">
        <link href="../css/fathmm.css" rel="stylesheet" type="text/css">
        
        <script type="text/javascript" src="../js/bootstrap.min.js"></script>
        <script type="text/javascript" src="../js/jquery.min.js"></script>
        
        <!-- AJAX Predictions -->
        <script type="text/javascript">
            var Refresh;
            var Container;

            $(document).ready(function()
            {
                Container = $("#Predictions");
                Refresh   = setInterval(
                    function()
                    {
                        Container.load('../tmp/""" + HTMLSession + """.htm');
                    }, 5000);
            });
        </script>
    </head>

    <body>
        <div class="container">
            <div class="navbar navbar-fixed-top">
                <div class="navbar-inner">
                    <div class="container">
                        <a class="btn btn-navbar" data-toggle="collapse" data-target=".nav-collapse">
                            <span class="icon-bar"></span>
                            <span class="icon-bar"></span>
                            <span class="icon-bar"></span>
                        </a>
                        <a class="brand" href="#">fathmm</a>
                        <div class="nav-collapse collapse">
                            <ul class="nav">
                                <li><a href="../index.html">Home</a></li>
                                <li><a href="https://github.com/HAShihab/fathmm">Download</a></li>
                                <li class="dropdown">
                                    <a href="#" class="dropdown-toggle" data-toggle="dropdown">New Submission <b class="caret"></b></a>
                                    <ul class="dropdown-menu">
                                        <li><a href="../inherited.html">Analyze dbSNP/Protein Variants</a></li>
                                        <li><a href="../cancer.html">Analyze Cancer-Associated Variants</a></li>
                                    </ul>
                                </li>
                            </ul>
                            <form class="navbar-form pull-right" action="./results.cgi" method="post" enctype="multipart/form-data">
                                <input class="span3" type="text" id="session" name="session" placeholder="Enter Your Job/Session Identifier">
                                <button type="submit" class="btn">Fetch Results</button>
                            </form>
                        </div><!--/.nav-collapse -->
                    </div>
                </div>
            </div>
            
            <div class="hero-unit">
                <h2>fathmm Predictions</h2>
                <p>
                    your request (Job/Session ID """ + HTMLSession + """)  is now being processed - if you submitted a large batch; please be patient as this may take a while.
                </p>
                <p>
                    <a href="#" class="btn btn-warning disabled btn-large pull-right" id="info" name="info">Processing Request ...</a>
                </p>
            </div>
            
            <div class="row">
                <div class="span12">
                    <div id="Predictions"></div>
                </div>
            </div>
                    
            <hr>
            <footer>
            <p>
            We welcome any comments and/or suggestions that you may have regarding our software and server - please send an email directly to fathmm@biocompute.org.uk
            </p>
            </footer>
        </div>
    
        <script type="text/javascript" src="../js/bootstrap.min.js"></script>
        <script type="text/javascript" src="../js/jquery.min.js"></script>
    </body>
</html>
        """
    #
    except Exception, e:
        print "Content-Type: text/html"
        print
        print e
