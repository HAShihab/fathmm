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
    
    HTML            = ""
    try:
        HTMLForm    = cgi.FieldStorage()
        HTMLSession = str(HTMLForm['session'].value).strip()
        
        HTML+= """
    <!DOCTYPE html>
        <head>
            <title>fathmm - fathmm Predictions</title>

            <meta charset="utf-8">
            
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <meta name="description" content="Predict the Functional and Phenotypic Consequences of Protein Missense Variants">
            <meta name="author" content="Hashem A. Shihab">
            
            <meta http-equiv="cache-control" content="max-age=0" />
            <meta http-equiv="cache-control" content="no-cache" />
            <meta http-equiv="expires" content="0" />
            <meta http-equiv="expires" content="Tue, 01 Jan 1980 1:00:00 GMT" />
            <meta http-equiv="pragma" content="no-cache" />
            
            <link rel="icon" type="image/png" href="../img/favicon.ico">
            
            <!-- HTML5 "Upscaling" -->
            <!--[if lt IE 9]>
                <script src="http://html5shiv.googlecode.com/svn/trunk/html5.js"></script>
            <![endif]-->

            <link href="../css/bootstrap.css" rel="stylesheet">
            <link href="../css/bootstrap-responsive.css" rel="stylesheet">
            <link href="../css/fathmm.css" rel="stylesheet" type="text/css">
            
            <script type="text/javascript" src="../js/bootstrap.min.js"></script>
            <script type="text/javascript" src="../js/jquery.min.js"></script>
            
            <!-- DATA TABLE -->
            <link href="../css/Data-Table.css" rel="stylesheet" type="text/css">
            
            <script type="text/javascript" src="../js/jquery.dataTables.min.js"></script>
            <script type="text/javascript" src="../js/Data-Table.js"></script>
            <!-- DATA TABLE -->
            
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
                                    <li><a href="../about.html">About</a></li>
                                    <li><a href="https://github.com/HAShihab/fathmm">Download</a></li>
                                </ul>
                                <form class="navbar-form pull-right" action="./results.cgi" method="post" enctype="multipart/form-data">
                                    <input class="span3" type="text" id="session" name="session" placeholder="Enter Your Job/Session Identifier">
                                    <button type="submit" class="btn">Fetch Results</button>
                                </form>
                            </div><!--/.nav-collapse -->
                        </div>
                    </div>
                </div>
        """
        
        
        if os.path.exists("../tmp/" + HTMLSession + ".txt"):
            # submission has been made ...
            if os.path.exists("../tmp/" + HTMLSession + ".tab"):
                # prediction have been made ...
                Predictions = [ x for x in open("../tmp/" + HTMLSession + ".tab", "r").readlines() if not x.startswith("#") ]
                
                HTML+= """
                <div class="hero-unit">
                    <h2>fathmm Predictions</h2>
                    <br />
                    <p>
                        Your request has now been processed and our predictions have been written onto your screen below - you can retrieve your results at 
                        any time during the next week using the following Job/Session Identifier: <strong>""" + HTMLSession + """</strong>. In addition, a 
                        "tab-delimited" report of our predictions has been made available via the "Download" button below.
                    </p>
                    <p>
                        <a href="../tmp/""" + HTMLSession + """.tab" class="btn btn-primary btn-large pull-right" id="info" name="info">Download &raquo;</a>
                    </p>
                </div>
                """
                
                if len(Predictions) > 1000:
                    HTML+= """
                <div class="row">
                    <div class="span10 offset1">
                        <div class="alert">
                            <button type="button" class="close" data-dismiss="alert">&times;</button>
                            <font color="#555555">
                                It appears that your submission exceeded 1,000 records; therefore, in order to minimize the load on your browser, we have limited the 
                                below to the first 1,000 records.  However, <strong>you can retrieve all your results in a tab-delimited report via the "Download" 
                                button above</strong> (right-click/save-as).
                            </font>
                        </div>
                    </div>
                </div>
                <br />
                    """
                HTML+= """
                <div class="row">
                    <div class="span12">
                        <table class="table table-striped table-bordered" id="results">
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
                
                for record in Predictions[:1000]:
                    if record and not record.startswith("#"):
                        record = record.split("\t")
                        
                        if record[4] in [ "DAMAGING", "CANCER" ]:
                            record[4] = "<p class='text-error'>%s</p>"% record[4]
                        else:
                            record[4] = "<p class='text-success'>%s</p>"% record[4]
                        
                        # warning/phenotypes formatting
                        if record[7].strip():
                            record[6]     = record[7]
                        else:
                            if record[6].strip():
                                record[6] = "<ul>" + "".join([ "<li>" + x + "</li>" for x in record[6].split("|") ]) + "</ul>"
                                
                        HTML+= """
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
                        
                HTML+= """
                            </tbody>
                        </table>
                    </div>
                </div>
                """
            else:
                # request is being processed ...
                HTML+= """
                <div class="hero-unit">
                    <h2>fathmm Predictions</h2>
                    <br />
                    <p>
                        Your request is now being processed and our predictions will be written onto your screen once completed - if you submitted a large batch, 
                        please be patient as this may take a while. Note, you can retrieve your results at any time during the next week using the following 
                        Job/Session Identifier: <strong>""" + HTMLSession + """</strong>. 
                    </p>
                    <p>
                        <a href="#" class="btn btn-warning disabled btn-large pull-right" id="info" name="info">Processing Request ...</a>
                    </p>
                </div>
                
                <script type="text/javascript">
                  setTimeout(function(){
                    location = ''
                  },2500)
                </script>
                """ 
        else:
            # no submission has been made/previous submission has been deleted ...
            HTML+= """
            <div class="hero-unit">
                <h2>fathmm Predictions</h2>
                <br />
                <p>
                    You have entered an invalid or unknown Job/Session Identifier - please note, predictions are stored on our server for one week before being deleted.
                </p>
            </div>
            """ 
        
        
        HTML+= """
            <hr>
            <footer>
            <p>
                If you have found this resource useful, please cite the following work:
                <br />
                <a href="http://www.ncbi.nlm.nih.gov/pubmed/23033316">Shihab HA, Gough J, Cooper DN, Stenson PD, Barker GLA, Edwards KJ, Day INM, Gaunt, TR. (2013). Predicting the Functional, Molecular and Phenotypic Consequences of Amino Acid Substitutions using Hidden Markov Models. Hum. Mutat., <b>34</b>, 57-65 </a>
            </p>
            <p>
                We welcome any comments and/or suggestions that you may have regarding our software and server - please send an email directly to fathmm@biocompute.org.uk
            </p>
        </footer>
        </div>
        
        <script type="text/javascript" src="../js/bootstrap.min.js"></script>
    </body>
</html>
    """
    except Exception, e:
        HTML = """
<!DOCTYPE html>
<!--[if lt IE 7]> <html class="lt-ie9 lt-ie8 lt-ie7"> <![endif]-->
<!--[if IE 7]>    <html class="lt-ie9 lt-ie8"> <![endif]-->
<!--[if IE 8]>    <html class="lt-ie9"> <![endif]-->
<!--[if gt IE 8]><!--> <html class=""> <!--<![endif]-->
    <head>
        <title>fathmm - fathmm Predictions</title>

        <meta charset="utf-8">
        
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <meta name="description" content="Predict the Functional and Phenotypic Consequences of Protein Missense Variants">
        <meta name="author" content="Hashem A. Shihab">
        
        <meta http-equiv="cache-control" content="max-age=0" /><script type="text/javascript" src="../js/jquery.dataTables.min.js"></script>
            <script type="text/javascript" src="../js/Data-Table.js"></script>
        <meta http-equiv="cache-control" content="no-cache" />
        <meta http-equiv="expires" content="0" />
        <meta http-equiv="expires" content="Tue, 01 Jan 1980 1:00:00 GMT" />
        <meta http-equiv="pragma" content="no-cache" />
        
        <link rel="icon" type="image/png" href="../img/favicon.ico">
        
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
                                <li><a href="../about.html">About</a></li>
                                <li><a href="https://github.com/HAShihab/fathmm">Download</a></li>
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
                <br />
                <p>
                    It appears there is a problem with our server - please try again and/or report the problem to the system administrator
                    using the following address: fathmm.biocompute.org.uk
                </p>
            </div>
            
            <div class="row">
                <div class="span12">
                    <form class="form-horizontal" id="myForm" name="myForm"">
                        <legend><h3>Stack Trace:</h3></legend>

                        <div class="row-fluid">
                            <div class="span12">
                            
                                <div class="well well-large">
                                    <div class="row-fluid">
                                        <div class="span12">
            """
        HTML+= traceback.format_exc()
        HTML+= """
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
                    If you have found this resource useful, please cite the following work:
                    <br />
                    <a href="http://www.ncbi.nlm.nih.gov/pubmed/23033316">Shihab HA, Gough J, Cooper DN, Stenson PD, Barker GLA, Edwards KJ, Day INM, Gaunt, TR. (2013). Predicting the Functional, Molecular and Phenotypic Consequences of Amino Acid Substitutions using Hidden Markov Models. Hum. Mutat., <b>34</b>, 57-65 </a>
                </p>
                <p>
                    We welcome any comments and/or suggestions that you may have regarding our software and server - please send an email directly to fathmm@biocompute.org.uk
                </p>
            </footer>
        </div>
        
        <script type="text/javascript" src="../js/bootstrap.min.js"></script>
    </body>
</html>
        """
    finally:
        print HTML
