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
        # if new submission, perform our analysis/predictions in the background ...
        HTMLForm    = cgi.FieldStorage()
        HTMLTime    = time.gmtime()
        HTMLSession = None
        
        if not HTMLForm.has_key('session'):
            while True:
                # assert a unique session is generated ...
                HTMLSession = str(uuid.uuid4())
                
                if not os.path.exists("../tmp/" + HTMLSession + ".txt"):
                    break
            
            # write the submission (& server parameters) to the server
            for x in [ "# IP: " + str(os.environ['REMOTE_ADDR']) + "\n",
                       "# Date: " + time.strftime("%d-%m-%Y", HTMLTime) + "\n",
                       "# Time: " + time.strftime("%H:%M:%S", HTMLTime) + "\n",
                       "# Weights: " + str(HTMLForm['weighted'].value).upper() + "\n",
                       "# Phenotypes: " + str(HTMLForm['phenotypes'].value).upper() + "\n",
                       str(HTMLForm['batch'].value)
                     ]:
                open("../tmp/" + HTMLSession + ".txt", "a").write(x)

            pid = subprocess.Popen([ sys.executable,
                                     'fathmm.py',
                                     '-i', '../tmp/' + HTMLSession + '.txt',
                                     '-o', '../tmp/' + HTMLSession + '.tab',
                                     '-w', str(HTMLForm['weighted'].value),
                                     '-p', str(HTMLForm['phenotypes'].value)
                                   ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        #
        else:
            # ...otherwise, fetch predictions from a previous session
            HTMLSession = str(HTMLForm['session'].value).strip()
        
        # assert the submission exists ...
        if os.path.exists("../tmp/" + HTMLSession + ".txt"):
            # assert results/predictions have been written
            if os.path.exists("../tmp/" + HTMLSession + ".tab"):
                print """
    <!DOCTYPE html>
    <!--[if lt IE 7]> <html class="lt-ie9 lt-ie8 lt-ie7"> <![endif]-->
    <!--[if IE 7]> <html class="lt-ie9 lt-ie8"> <![endif]-->
    <!--[if IE 8]> <html class="lt-ie9"> <![endif]-->
    <!--[if gt IE 8]><!--> <html class=""> <!--<![endif]-->
        <head>
            <title>fathmm - fathmm Predictions</title>

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
            
            <script type="text/javascript">
              var _gaq = _gaq || [];
              _gaq.push(['_setAccount', 'UA-29568329-3']);
              _gaq.push(['_trackPageview']);

              (function() {
                var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
                ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
                var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
              })();
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
                            <a class="brand" href="">fathmm</a>
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
                        Your request has now been processed and our predictions have been written onto your screen - a "tab-delimited" report of our predictions
                        has been made available below.&nbsp;&nbsp;In addition, you can retrieve your results at any time using the below Job/Session Identifier.
                        <br />
                        <br />
                        Job/Session ID: """ + HTMLSession + """
                    </p>
                    <p>
                        <a href="../tmp/""" + HTMLSession + """.tab" class="btn btn-primary btn-large pull-right" id="info" name="info">Download &raquo;</a>
                    </p>
                </div>
                
                <div class="row">
                    <div class="span12">
                        <table class="table table-striped">
                            <thead>
                                <tr>
                                    <th>#</th>
                                    <th>dbSNP ID</th>
                                    <th>Protein ID</th>
                                    <th>Substitution</th>
                                    <th>Prediction</th>
                                    <th>Score</th>
                                    <th>Further Information</th>
                                </tr>
                            </thead>
                            <tbody>
                """
                
                for record in open("../tmp/" + HTMLSession + ".tab", "r"):
                    if record and not record.startswith("#"):
                        record = record.split("\t")
                        
                        # prediction formatting
                        if record[4] in [ "DAMAGING", "CANCER" ]:
                            record[4] = "<p class='text-error'>%s</p>"% record[4]
                        else:
                            record[4] = "<p class='text-success'>%s</p>"% record[4]
                        
                        # warning/phenotypes formatting
                        if record[7].strip():
                            record[6] = record[7]
                        else:
                            if record[6].strip():
                                record[6] = "<ul>" + "".join([ "<li>" + x + "</li>" for x in record[6].split("|") ]) + "</ul>"
                                
                        print """
                                <tr>
                                    <td>""" + record[0] + """</td>
                                    <td>""" + record[1] + """</td>
                                    <td>""" + record[2] + """</td>
                                    <td>""" + record[3] + """</td>
                                    <td>""" + record[4] + """</td>
                                    <td>""" + record[5] + """</td>
                                    <td>""" + record[6] + """</td>
                                </tr>
                        """
                        
                print """
                            </tbody>
                        </table>
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
            else:
                print """
<!DOCTYPE html>
<!--[if lt IE 7]> <html class="lt-ie9 lt-ie8 lt-ie7"> <![endif]-->
<!--[if IE 7]>    <html class="lt-ie9 lt-ie8"> <![endif]-->
<!--[if IE 8]>    <html class="lt-ie9"> <![endif]-->
<!--[if gt IE 8]><!--> <html class=""> <!--<![endif]-->
    <head>
        <title>fathmm - Processing Submission</title>

        <meta charset="utf-8">
        
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <meta name="description" content="Predict the Functional and Phenotypic Consequences of Protein Missense Variants">
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
        
        <script type="text/JavaScript">
            function redirect() {
                window.location = './results.cgi?session=""" + HTMLSession + """';
            }
        </script>
        
        <script type="text/javascript">
          var _gaq = _gaq || [];
          _gaq.push(['_setAccount', 'UA-29568329-3']);
          _gaq.push(['_trackPageview']);

          (function() {
            var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
            ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
            var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
          })();
        </script>
    </head>

    <body onLoad="setTimeout('redirect()', 5000)">
        <div class="container">
        
            <div class="navbar navbar-fixed-top">
                <div class="navbar-inner">
                    <div class="container">
                        <a class="btn btn-navbar" data-toggle="collapse" data-target=".nav-collapse">
                            <span class="icon-bar"></span>
                            <span class="icon-bar"></span>
                            <span class="icon-bar"></span>
                        </a>
                        <a class="brand" href="">fathmm</a>
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
                    Your request is now being processed and our predictions will be written onto your screen once completed - if you submitted a large batch, 
                    please be patient as this may take a while.&nbsp;&nbsp;Alternatively, you can retrieve your results at any time using the below Job/Session 
                    Identifier.
                    <br />
                    <br />
                    Job/Session ID: """ + HTMLSession + """
                </p>
                <p>
                    <a href="#" class="btn btn-warning disabled btn-large pull-right" id="info" name="info">Processing Request ...</a>
                </p>
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
        else:
            # ... either there is no submission or the previous submission 
            # has been deleted
            print """
<!DOCTYPE html>
<!--[if lt IE 7]> <html class="lt-ie9 lt-ie8 lt-ie7"> <![endif]-->
<!--[if IE 7]>    <html class="lt-ie9 lt-ie8"> <![endif]-->
<!--[if IE 8]>    <html class="lt-ie9"> <![endif]-->
<!--[if gt IE 8]><!--> <html class=""> <!--<![endif]-->
    <head>
        <title>fathmm - Unknown Job/Session Identifier</title>

        <meta charset="utf-8">
        
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <meta name="description" content="Predict the Functional and Phenotypic Consequences of Protein Missense Variants">
        <meta name="author" content="Hashem A. Shihab">
        
        <meta http-equiv="cache-control" content="max-age=0" />
        <meta http-equiv="cache-control" content="no-cache" />
        <meta http-equiv="expires" content="0" />
        <meta http-equiv="expires" content="Tue, 01 Jan 1980 1:00:00 GMT" />
        <meta http-equiv="pragma" content="no-cache" />
        
        <!--[if lt IE 9]>
            <script src="http://html5shiv.googlecode.com/svn/trunk/html5.js"></script>
        <![endif]-->

        <link href="../css/bootstrap.css" rel="stylesheet">
        <link href="../css/bootstrap-responsive.css" rel="stylesheet">
        <link href="../css/fathmm.css" rel="stylesheet" type="text/css">
        
        <script type="text/javascript" src="../js/bootstrap.min.js"></script>
        <script type="text/javascript" src="../js/jquery.min.js"></script>
        
        <script type="text/javascript">
          var _gaq = _gaq || [];
          _gaq.push(['_setAccount', 'UA-29568329-3']);
          _gaq.push(['_trackPageview']);

          (function() {
            var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
            ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
            var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
          })();
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
                        <a class="brand" href="">fathmm</a>
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
                <h2>Unknown Job/Session Identifier</h2>
                <p>
                    You have entered an invalid or unknown Job/Session Identifier - please note, predictions are typically stored on our server for one week before being deleted.
                </p>
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
        print """
<!DOCTYPE html>
<!--[if lt IE 7]> <html class="lt-ie9 lt-ie8 lt-ie7"> <![endif]-->
<!--[if IE 7]>    <html class="lt-ie9 lt-ie8"> <![endif]-->
<!--[if IE 8]>    <html class="lt-ie9"> <![endif]-->
<!--[if gt IE 8]><!--> <html class=""> <!--<![endif]-->
    <head>
        <title>fathmm - Server Error</title>

        <meta charset="utf-8">
        
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <meta name="description" content="Predict the Functional and Phenotypic Consequences of Protein Missense Variants">
        <meta name="author" content="Hashem A. Shihab">
        
        <meta http-equiv="cache-control" content="max-age=0" />
        <meta http-equiv="cache-control" content="no-cache" />
        <meta http-equiv="expires" content="0" />
        <meta http-equiv="expires" content="Tue, 01 Jan 1980 1:00:00 GMT" />
        <meta http-equiv="pragma" content="no-cache" />
        
        <!--[if lt IE 9]>
            <script src="http://html5shiv.googlecode.com/svn/trunk/html5.js"></script>
        <![endif]-->

        <link href="../css/bootstrap.css" rel="stylesheet">
        <link href="../css/bootstrap-responsive.css" rel="stylesheet">
        <link href="../css/fathmm.css" rel="stylesheet" type="text/css">
        
        <script type="text/javascript" src="../js/bootstrap.min.js"></script>
        <script type="text/javascript" src="../js/jquery.min.js"></script>
        
        <script type="text/javascript">
          var _gaq = _gaq || [];
          _gaq.push(['_setAccount', 'UA-29568329-3']);
          _gaq.push(['_trackPageview']);

          (function() {
            var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
            ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
            var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
          })();
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
                        <a class="brand" href="">fathmm</a>
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
                <h2>Server Error</h2>
                <p>
                    It appears there is a problem with our server - please try again and/or report the problem to the system administrator
                    using the following address: fathmm.biocompute.org.uk
                </p>
            </div>
            
            <div class="row">
                <div class="span12">
                    <form class="form-horizontal" id="myForm" name="myForm" action="./cgi-bin/results.cgi" method="post" enctype="multipart/form-data">
                        <legend><h3>Stack Trace:</h3></legend>

                        <div class="row-fluid">
                            <div class="span12">
                            
                                <div class="well well-large">
                                    <div class="row-fluid">
                                        <div class="span12">
            """
        
        print traceback.format_exc()

        print """
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>      
                    </form>
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
