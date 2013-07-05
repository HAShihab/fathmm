#!/usr/bin/python -u

import re
import math
import argparse
import ConfigParser

import MySQLdb
import MySQLdb.cursors

#
def map_position(domain, substitution):
    """
    """
    
    if int(substitution[1:-1]) < int(domain['seq_begin']) or \
       int(substitution[1:-1]) > int(domain['seq_end']):
           return None
    
    x = int(domain['seq_begin']) - 1 
    y = int(domain['hmm_begin']) - 1 
    
    for residue in list(domain['align']):
        if residue.isupper() or residue.islower():
            # upper-case (match)/lower-case (insertion) characters correspond 
            # to a sequence residue
            x += 1
        if residue.isupper() or residue == "-":
            # upper-case (match)/gap (deletion) characters correspond to
            # a HMM position
            y += 1
        if x == int(substitution[1:-1]):
            if residue.isupper():
                return str(y)
            # unmapped residue (i.e. insertion)
            return None
    return None

#
def fetch_phenotype_prediction(Facade, Phenotype):
    """
    """
    
    Phenotypes = ""
    
    if not Arg.phenotypes:
        return Phenotypes
    
    for x in sorted(Facade, key=lambda x:x['information'], reverse=True):
        # use the most informative HMM
        if x['accession']:
            dbCursor.execute("select * from PHENOTYPES where accession='" + x['accession'] + "' and type='" + Phenotype + "' and origin=1 order by score")
            Phenotypes = "|".join([ x['description'] for x in dbCursor.fetchall() ])
            
            if Phenotypes:
                break
    
    return Phenotypes

#
def Process(Protein, Substitution, Weights=None, Cutoff=None, Phenotype=None):
    """
    """
    
    Processed = {
                 'Prediction' : "",
                 'Score':       "",       
                 'Phenotype':   "",
                 'HMM':         "",
                 'Description': "",
                 'Position':    "",
                 'W':           "",
                 'M':           "",
                 'D':           "",
                 'O':           "",
                 'Warning':     "",
                }
    
    # fetch pre-computed sequence record
    dbCursor.execute("select a.* from SEQUENCE a, PROTEIN b where a.id=b.id and b.name='" + Protein + "'")
    Sequence = dbCursor.fetchone()
    
    if not Sequence:
        Processed['Warning'] = "No Sequence Record Found"
        
        return Processed
    
    # authenticate protein/substitution
    Warning     = None
    if not Warning and not re.compile("^[ARNDCEQGHILKMFPSTWYV]\d+[ARNDCEQGHILKMFPSTWYV]$", re.IGNORECASE).match(Substitution):
        Warning = "Invalid Substitution Format"
    if not Warning and int(Substitution[1:-1]) > len(Sequence['sequence']):
        Warning = "Invalid Substitution Position"
    if not Warning and not Substitution[0] == Sequence['sequence'][int(Substitution[1:-1]) - 1]:
        Warning = "Inconsistent Wild-Type Residue (Expected '" + Sequence['sequence'][int(Substitution[1:-1]) - 1] + "')"
    if not Warning and Substitution[0] == Substitution[-1]:
        Warning = "Synonymous Mutation"
        
    if Warning:
        Processed['Warning'] = Warning
        
        return Processed
    
    # fetch pre-computed domain assignment(s)
    dbCursor.execute("select * from DOMAINS where id='" + str(Sequence['id']) + "' and " + Substitution[1:-1] + " between seq_begin and seq_end order by score")
    Domain = dbCursor.fetchall()
    
    Facade = []
    
    for x in Domain:
        residue = map_position(x, Substitution)
        
        if residue:
            dbCursor.execute("select a.*, b.* from PROBABILITIES a, LIBRARY b where a.id=b.id and a.id='" + str(x['hmm']) + "' and a.position='" + residue + "'")
            Prob = dbCursor.fetchone()
            
            if Prob:
                Facade.append(Prob)
    
    #
    
    # fetch phenotype association(s)
    if Phenotype:
        Processed['Phenotype'] = fetch_phenotype_prediction(Facade, Phenotype)
    
    # derive/return a prediction ...
    if not Weights or Weights == "UNWEIGHTED":
        # append non-domain-based/sequence conservation to the Facade
        dbCursor.execute("select a.*, b.*  from PROBABILITIES a, LIBRARY b where a.id=b.id and a.id='" + str(Sequence['id']) + "' and a.position='" + Substitution[1:-1] + "'")
        Prob = dbCursor.fetchone()
        
        if Prob:
            Facade.append(Prob)
        
        for x in sorted(Facade, key=lambda x:x['information'], reverse=True):
            try:
                Processed['HMM']         = x['id']
                Processed['Description'] = x['description']
                Processed['Position']    = x['position']
                Processed['W']           = x[Substitution[0]]
                Processed['M']           = x[Substitution[-1]]
                Processed['D']           = ""
                Processed['O']           = ""
                
                Processed['Score']       = "%.2f" % math.log((Processed['M'] / (1.0 - Processed['M'])) / (Processed['W'] / (1.0 - Processed['W'])), 2)

                Processed['Prediction']  = ""                
                if Cutoff:
                    if float(Processed['Score']) <= Cutoff: Processed['Prediction'] = "DAMAGING"
                    if float(Processed['Score']) >  Cutoff: Processed['Prediction'] = "TOLERATED"
                
                return Processed
            except Exception, e:
                # skip the rare occasion(s) when probabilities are zero
                pass
    else:
        # prioritize domain-based predictions
        for x in sorted(Facade, key=lambda x:x['information'], reverse=True):
            try:
                dbCursor.execute("select * from WEIGHTS where id='" + x['id'] + "' and type='" + Weights + "'")
                w = \
                    dbCursor.fetchone()
                
                if w:
                    Processed['HMM']         = x['id']
                    Processed['Description'] = x['description']
                    Processed['Position']    = x['position']
                    Processed['W']           = x[Substitution[0]]
                    Processed['M']           = x[Substitution[-1]]
                    Processed['D']           = w['disease'] + 1.0
                    Processed['O']           = w['other'] + 1.0
                    
                    Processed['Score']       = "%.2f" % math.log(((1.0 - Processed['W']) * Processed['O']) / ((1.0 - Processed['M']) * Processed['D']), 2)

                    Processed['Prediction']  = ""                
                    if Cutoff:
                        if Weights == "INHERITED":
                            if float(Processed['Score']) <= Cutoff: Processed['Prediction'] = "DAMAGING"
                            if float(Processed['Score']) >  Cutoff: Processed['Prediction'] = "TOLERATED"
                        
                        if Weights == "CANCER":
                            if float(Processed['Score']) <= Cutoff: Processed['Prediction'] = "CANCER"
                            if float(Processed['Score']) >  Cutoff: Processed['Prediction'] = "PASSENGER/OTHER"
                    
                    return Processed
            except Exception, e:
                # skip the rare occasion(s) when probabilities are zero
                pass
        
        # no domain-based prediction, use a non-domain-based/sequence conservation prediction
        dbCursor.execute("select a.*, b.*  from PROBABILITIES a, LIBRARY b where a.id=b.id and a.id='" + str(Sequence['id']) + "' and a.position='" + Substitution[1:-1] + "'")
        Facade = dbCursor.fetchone()
        
        if Facade:
            try:
                dbCursor.execute("select * from WEIGHTS where id='" + Facade['id'] + "' and type='" + Weights + "'")
                w = \
                    dbCursor.fetchone()
                
                if w:
                    Processed['HMM']         = Facade['id']
                    Processed['Description'] = Facade['description']
                    Processed['Position']    = Facade['position']
                    Processed['W']           = Facade[Substitution[0]]
                    Processed['M']           = Facade[Substitution[-1]]
                    Processed['D']           = w['disease'] + 1.0
                    Processed['O']           = w['other'] + 1.0
                    
                    Processed['Score']       = "%.2f" % math.log(((1.0 - Processed['W']) * Processed['O']) / ((1.0 - Processed['M']) * Processed['D']), 2)

                    Processed['Prediction']  = ""                
                    if Cutoff:
                        if Weights == "INHERITED":
                            if float(Processed['Score']) <= Cutoff: Processed['Prediction'] = "DAMAGING"
                            if float(Processed['Score']) >  Cutoff: Processed['Prediction'] = "TOLERATED"
                        
                        if Weights == "CANCER":
                            if float(Processed['Score']) <= Cutoff: Processed['Prediction'] = "CANCER"
                            if float(Processed['Score']) >  Cutoff: Processed['Prediction'] = "PASSENGER/OTHER"
                    
                    return Processed
            except Exception, e:
                # skip the rare occasion(s) when probabilities are zero
                raise        

    return None



#
if __name__ == '__main__':
    #
    
    Config   = ConfigParser.ConfigParser()
    Config.read("./config.ini")
        
    dbCursor = MySQLdb.connect(
                   host     = str(Config.get("DATABASE", "HOST")),
                   port     = int(Config.get("DATABASE", "PORT")),
                   user     = str(Config.get("DATABASE", "USER")),
                   passwd   = str(Config.get("DATABASE", "PASSWD")),
                   db       = str(Config.get("DATABASE", "DB")),
                   compress = 1
               ).cursor(MySQLdb.cursors.DictCursor)
           
    # parse program argument(s)
    parser = argparse.ArgumentParser(
                                     description = 'Functional Analysis through Hidden Markov Models',
                                     add_help    = False
                                    )
    parser.add_argument(
                        "-h",
                        "--help",
                        action = "help",
                        help   = argparse.SUPPRESS
                       )
    
    group = \
        parser.add_argument_group("Required")
    group.add_argument(
                       'fi',
                       metavar = '<F1>', 
                       type    = argparse.FileType("r"),
                       help    = 'a file containing the mutation data to process'
                      )
    group.add_argument(
                       'fo',
                       metavar = '<F2>',
                       type    = argparse.FileType("w"),
                       help    = 'where predictions/phenotype-associations will be written'
                      )
    
    group = \
        parser.add_argument_group("Options")
    group.add_argument(
                       '-w',
                       dest    = 'weights',
                       metavar = "<S>",
                       default = "INHERITED",
                       help    = "use pathogenicity weights <S> when returning predictions"
                      )
    group.add_argument(
                       '-t',
                       dest    = 'threshold',
                       metavar = "<N>",
                       default = None,
                       type    = float,
                       help    = "use prediction threshold <N> when returning predictions"
                      )
    group.add_argument(
                       '-p',
                       dest    = 'phenotypes',
                       metavar = "<S>",
                       default = "",
                       help    = "use phenotype ontology <S> when returning domain-phenotype associations"
                      ); Arg = parser.parse_args()
    
    if Arg.weights:
        # authenticate prediction weights
        dbCursor.execute("select distinct type from WEIGHTS")
        
        if not Arg.weights.upper() in [ x['type'] for x in dbCursor.fetchall() ] + [ "UNWEIGHTED" ]: 
            parser.error("argument -w: invalid option: '{0}'".format(Arg.weights))
        
        if Arg.threshold == None:
            # initialize predcition threshold
            if Arg.weights.upper() == "UNWEIGHTED": Arg.threshold = -3.00
            if Arg.weights.upper() == "INHERITED":  Arg.threshold = -1.50
            if Arg.weights.upper() == "CANCER":     Arg.threshold = -0.75
    
    if Arg.phenotypes:
        # authenticate prediction phenotype
        dbCursor.execute("select distinct type from PHENOTYPES")
        
        if not Arg.phenotypes.upper() in [ x['type'] for x in dbCursor.fetchall() ]:
            parser.error("argument -p: invalid option: '{0}'".format(Arg.phenotypes))

    
    ##


    Arg.fo.write("\t".join([ "#", 
                             "dbSNP ID",
                             "Protein ID",
                             "Substitution",
                             "Prediction",
                             "Score",
                             "Domain-Phenotype Association",
                             "Warning",
                             "HMM ID",
                             "HMM Description",
                             "HMM Pos.",
                             "HMM Prob. W.",
                             "HMM Prob. M.",
                             "HMM Weights D.",
                             "HMM Weights O."
                           ]) + "\n")
    
    idx = 1
    for record in Arg.fi:
        record = record.strip()
        
        if record and not record.startswith("#"):
            try:
                if re.compile("^rs\d+$", re.IGNORECASE).match(record):
                    dbCursor.execute("select distinct * from VARIANTS where id='" + record + "'")
                    dbRecords = dbCursor.fetchall()
                    
                    if not dbRecords:
                        Arg.fo.write(
                            "\t".join([ str(idx),
                                        record,
                                        "",
                                        "",
                                        "",
                                        "",
                                        "",
                                        "No dbSNP Mapping(s)",
                                        "",
                                        "",
                                        "",
                                        "",
                                        "",
                                        "",
                                        ""
                                      ]) + "\n"
                        ); idx += 1; continue
                        
                    for x in dbRecords:
                        dbSNP        = x['id']
                        Protein      = x['protein']
                        Substitution = x['substitution']
                        
                        Prediction = \
                            Process(Protein, Substitution, Weights=Arg.weights.upper(), Cutoff=Arg.threshold, Phenotype=Arg.phenotypes.upper())

                        if not Prediction:
                            Arg.fo.write(
                                "\t".join([ str(idx),
                                            dbSNP,
                                            Protein,
                                            Substitution,
                                            "",
                                            "",
                                            "",
                                            "No Prediction Available",
                                            "",
                                            "",
                                            "",
                                            "",
                                            "",
                                            "",
                                            ""
                                          ]) + "\n"
                            ); idx += 1; continue
                        Arg.fo.write(
                            "\t".join([ str(idx),
                                        dbSNP,
                                        Protein,
                                        Substitution,
                                        str(Prediction['Prediction']),
                                        str(Prediction['Score']),
                                        str(Prediction['Phenotype']),
                                        str(Prediction['Warning']),
                                        str(Prediction['HMM']),
                                        str(Prediction['Description']),
                                        str(Prediction['Position']),
                                        str(Prediction['W']),
                                        str(Prediction['M']),
                                        str(Prediction['D']),
                                        str(Prediction['O']) 
                                      ]) + "\n"
                        ); idx += 1; continue
                else:
                    dbSNP         = ""
                    Protein       = record.upper().split()[0]

                    for Substitution in [ x.strip() for x in record.upper().split()[1].split(",") ]:
                        Prediction = \
                            Process(Protein, Substitution, Weights=Arg.weights.upper(), Cutoff=Arg.threshold, Phenotype=Arg.phenotypes.upper())

                        if not Prediction:
                            Arg.fo.write(
                                "\t".join([ str(idx),
                                            dbSNP,
                                            Protein,
                                            Substitution,
                                            "",
                                            "",
                                            "",
                                            "No Prediction Available",
                                            "",
                                            "",
                                            "",
                                            "",
                                            "",
                                            "",
                                            ""
                                          ]) + "\n"
                            ); idx += 1; continue
                        Arg.fo.write(
                            "\t".join([ str(idx),
                                        dbSNP,
                                        Protein,
                                        Substitution,
                                        str(Prediction['Prediction']),
                                        str(Prediction['Score']),
                                        str(Prediction['Phenotype']),
                                        str(Prediction['Warning']),
                                        str(Prediction['HMM']),
                                        str(Prediction['Description']),
                                        str(Prediction['Position']),
                                        str(Prediction['W']),
                                        str(Prediction['M']),
                                        str(Prediction['D']),
                                        str(Prediction['O']) 
                                      ]) + "\n"
                        ); idx += 1; continue
            except Exception, e:
                Arg.fo.write(
                    "\t".join([ str(idx),
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "An Error Occured While Processing The Record: " + record,
                                "",
                                "",
                                "",
                                "",
                                "",
                                ""
                              ]) + "\n"
                ); idx += 1

