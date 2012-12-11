fathmm
======

This is the source code for the Functional Analysis through Hidden Markov Models (fathmm) software and server.

## General Requirements

You will need the following packages installed on your system:

* MySQL
* Python & Python MySQLdb (tested with Python 2.6/2.7)

## Setup:

### Database

* our pre-computed database, including instructions on how to create/upload the database, can be found at ftp://supfam2.cs.bris.ac.uk/FATHMM/database
* create a configuration file named "config.ini" and enter the following (substituting the required information with your credentials):

	[DATABASE]
	HOST    = [MySQL Host]
	PORT    = [MySQL Port]
	USER    = [MySQL Username]
	PASSWD  = [MySQL Password]
	DB      = fathmm

### Software

* download ./cgi-bin/fathmm.py (the --help parameter can be used to view expected program parameters)

## Known Issues:

We welcome any comments and/or suggestions that you may have regarding our software and server - please send an email directly to fathmm@biocompute.org.uk
