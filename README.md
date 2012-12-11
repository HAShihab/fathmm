fathmm
======

This is the source code for the Functional Analysis through Hidden Markov Models (fathmm) software and server.

## Installing our software:

1. download/install our precomputed database following the instructions on ftp://supfam2.cs.bris.ac.uk/FATHMM/database/
2. create the configuration file "config.ini" and write the following (substituting the required information with your credentials):

    [DATABASE]
    HOST    = [MySQL Host]
    PORT    = [MySQL Port]
    USER    = [MySQL Username]
    PASSWD  = [MySQL Password]
    DB      = fathmm

3. download ./cgi-bin/fathmm.py and use the --help parameter for further information

## Known Issues:

We welcome any comments and/or suggestions that you may have regarding our software and server - please send an email directly to fathmm@biocompute.org.uk
