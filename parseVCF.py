#!/usr/bin/python -u

import os
import argparse

#
if __name__ == '__main__':
    #
    
    parser = argparse.ArgumentParser(
        description = """
This program will parse the annotations provided by the Ensembl Variant Predictor (VEP)
and create a list of protein consequences which can then be submitted to the 
Functional Analysis through Hidden Markov Models (FATHMM) software and server.
                      """,
        epilog      = """
Note: in order for this script to work, the VEP must be called using the additional --protein parameter.
                      """,
                                    )
    
    group = parser.add_argument_group('required arguments')
    group.add_argument(
                      "-i",
                      dest     = "input",
                      help     = "a file containing the VEP annotations to parse", 
                      metavar  = "<INPUT>",
                      default  = None,
                      required = True
                      )
    group.add_argument(
                      "-o",
                      dest     = "output",
                      help     = "where the protein consequences should be written to", 
                      metavar  = "<OUTPUT>",
                      default  = None,
                      required = True
                      )
    
    arguments = parser.parse_args()
    
    #

    try:
        Consequence = {}
        
        for record in open(arguments.input, "r"):
            if record.startswith("#") or not record.find("missense_variant") >= 0:
                continue
            record  = record.strip().split("\t")
            
            x       = record[10].split("/")
            Var     = x[0] + record[9] + x[1]
            Info    = dict([ x.split("=") for x in record[-1].split(";") ])
            
            if not Info.has_key("ENSP"):
                continue
            
            if not Consequence.has_key(Info['ENSP']):
                Consequence[Info['ENSP']] = []
            Consequence[Info['ENSP']].append(Var)
        
        Processed   = open(arguments.output, "w")
        for id in Consequence:
            Processed.write(id + " " + ",".join(set(Consequence[id])) + "\n")
    #
    except:
        raise
