fathmm
======

This is the source code for the Functional Analysis through Hidden Markov Models 
(fathmm) software and server.

## General Requirements

You will need the following packages installed on your system:

* MySQL
* Python & Python MySQLdb (tested with Python 2.6/2.7)

## Setup:

* our pre-computed database, including instructions on how to create/upload the 
database, can be found at ftp://supfam2.cs.bris.ac.uk/FATHMM/database
* create a configuration file named "config.ini" and enter the following (substituting 
the required information with your credentials):

```
[DATABASE]
HOST   = [MySQL Host]
PORT   = [MySQL Port]
USER   = [MySQL Username]
PASSWD = [MySQL Password]
DB     = fathmm
```

* download "fathmm.py" from ./cgi-bin

## Running our Software

In it's simplest form, our software parses dbSNP rs IDs and protein missense 
mutations from <file> and returns a list of predictions weighted for 
inherited-disease mutations (Human).  Furthermore, we return predictions
on disease-associations when a mutation falls within a SUPERFAMILY domain.

```
python fathmm.py <FILE> <OPTIONS>
```

The --help parameter can be used to view expected program parameters.  In brief,
there are two optional parameters which control the prediction algorithm and
disease-associations reported:

* -w [?]

This parameter controls the prediction algorithm used in our software:

```
Inherited  : return predictions weighted for human inherited-disease mutations [this is the default]
Unweighted : return unweighted (species-independant) predictions 
Cancer     : return predictions weighted specifically for Human cancer
```

* -p [?]

This parameter controls which phenotype ontology to use:

```
DO : Disease Ontology
HP : Human Phenotype Ontology
MP : Mouse Phenotype Ontology
WP : Worm Phenotype Ontology
YP : Yeast Phenotype Ontology
FP : Fly Phenotype Ontology
FA : Fly Anatomy Ontology
ZA : Zebrafish Anatomy Ontology
AP : Arabidopsis Plant Ontology
KW : UniProtKB KeyWords
```

## Known Issues:

We welcome any comments and/or suggestions that you may have regarding our software and server - please send an email directly to fathmm@biocompute.org.uk
