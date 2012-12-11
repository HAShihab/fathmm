#!/usr/bin/python -u

import os
import re
import math
import ConfigParser

import MySQLdb
import MySQLdb.cursors

from optparse import OptionParser, SUPPRESS_HELP

#
def map_position(domain, substitution):
    # return if the amino acid substitution falls outside of the domain 
    # assignment 
    if int(substitution[1:-1]) < int(domain['seq_start']) or \
       int(substitution[1:-1]) > int(domain['seq_end']):
        return None
    
    # adjust sequence/HMM position (Python has zero-based coordinates)
    x = int(domain['seq_start']) - 1 
    y = int(domain['hmm_start']) - 1 
    
    # map amino acid substitution onto the HMM
    for residue in list(domain['alignment']):
        if residue.isupper() or residue.islower():
            # upper-case (match)/lower-case (insertion) characters map 
            # onto a sequence residue
            x += 1
        if residue.isupper() or residue == "-":
            # upper-case (match)/gap (deletion) characters map onto a
            # HMM position
            y += 1
        if x == int(substitution[1:-1]):
            if residue.isupper():
                # return corresponding position within the HMM
                return str(y)
            # unmapped residue (i.e. insertion)
            return None
    return None
    
#
def process_record(dbSNP, protein, substitution):
    
    #
    # FETCH PRE-COMPUTED SEQUENCE RECORD
    #
    
    dbCursor.execute("select a.* from Clusters a, Protein b where a.id=b.cluster and b.name='" + protein + "'")
    SeqRecord = dbCursor.fetchone()
    
    if not SeqRecord:
        # no pre-computed sequence record ...
        Predictions.write(
            "\t".join([ str(idx), dbSNP, protein, substitution, "", "", "", "No Sequence Record Found" ]) + "\n"
        ); return False
    
    #
    # AUTHENTICATE PROTEIN/SUBSTITUTION
    #
    
    _         = None
    
    if not _ and not re.compile("^[ARNDCEQGHILKMFPSTWYV]\d+[ARNDCEQGHILKMFPSTWYV]$", re.IGNORECASE).match(substitution):
        _     = "Invalid Substitution Format"
    if not _ and int(substitution[1:-1]) > SeqRecord['sequence'].__len__():
        _     = "Invalid Substitution Position"
    if not _ and not substitution[0] == SeqRecord['sequence'][int(substitution[1:-1]) - 1]:
        _     = "Inconsistent Wild-Type Residue (Expected '" + SeqRecord['sequence'][int(substitution[1:-1]) - 1] + "')"
    if not _ and substitution[0] == substitution[-1]:
        _     = "Synonymous Mutation"
    if  _:
        Predictions.write(
            "\t".join([ str(idx), dbSNP, protein, substitution, "", "", "", _ ]) + "\n"
        ); return False
    
    #
    # CREATE A FACADE OF POSSIBLE PREDICTION(S)
    #
    
    Facade    = []
    
    # fetch substitution-harbouring protein domains ...
    dbCursor.execute("select * from Domains where id='" + SeqRecord['id'] + "' and " + substitution[1:-1] + " between seq_start and seq_end order by evalue")
    DomRecord = dbCursor.fetchall()
    
    for x in DomRecord:
        # map the substitution onto the HMMs
        residue = map_position(x, substitution)
            
        if residue:
            # fetch description/probabilities for mapped position
            dbCursor.execute("select a.*, b.accession, b.description from Probabilities a, Library b where a.id=b.id and a.id='" + str(x['hmm']) + "' and a.position='" + residue + "'")
            dbRecord = dbCursor.fetchone()
            
            if dbRecord:
                Facade.append(dbRecord)
    
    # ... append JackHMMER probabilities
    dbCursor.execute("select a.*, b.accession, b.description from Probabilities a, Library b where a.id=b.id and a.id='" + SeqRecord['id'] + "' and a.position='" + substitution[1:-1] + "'")
    ConRecord = dbCursor.fetchone()
    if ConRecord:
        Facade.append(ConRecord)
    
    if not Facade:
        # no domain assignments/JackHMMER conservation
        Predictions.write(
            "\t".join([ str(idx), dbSNP, protein, substitution, "", "", "", "No Prediction Available" ]) + "\n"
        ); return False
    
    #
    # DERIVE A PREDICTION USING THE MOST INFORMATIVE HMM WITHIN THE FACADE
    #
    
    Score      = None
    Phenotypes = None
    
    for x in sorted(Facade, key=lambda x:x['information'], reverse=True):
        if not Score:
            if options.weights.upper() == "UNWEIGHTED":
                # "Unweighted Prediction" ...
                p     = x[substitution[0]]
                q     = x[substitution[-1]]
                
                Score = "%.2f" % math.log((q / (1.0 - q)) / (p / (1.0 - p)), 2)
            else:
                # "Inherited/Cancer Prediction" ...
                dbCursor.execute("select * from Weights where id='" + x['id'] + "' and type='" + options.weights + "'")
                Weights   = dbCursor.fetchone()

                if Weights:
                    p     = x[substitution[0]]
                    q     = x[substitution[-1]]
                    r     = Weights['polymorphic']
                    s     = Weights['pathogenic']
                    
                    Score = "%.2f" % math.log(((1.0 - p) * (r + 1.0)) / ((1.0 - q) * (s + 1.0)), 2)
        
        if not Phenotypes:
            # domain-phenotype prediction (via SUPERFAMILY)
            dbCursor.execute("select * from Phenotypes where accession='" + x['accession'] + "' and ontology='" + options.phenotypes + "' and origin=1 order by score")
            Phenotypes = dbCursor.fetchall()
            
        if Score and Phenotypes:
            # prediction(s) made ...
            break
    
    #
    # WRITE OUR PREDICTION
    #
    
    if not Score:
        # no pathogenicity weights/prediction score
        Predictions.write(
            "\t".join([ str(idx), dbSNP, protein, substitution, "", "", "", "No Prediction Available" ]) + "\n"
        ); return False
    
    Tag = None

    if not options.weights.upper() == "UNWEIGHTED":
        if options.weights.upper() == "INHERITED":
            # "Inherited Disease" predictions ...
            Tag = "TOLERATED"
            
            if float(Score) < -1.50:
                Tag = "DAMAGING" 
        else:
            # "Cancer-Associated" predictions ...
            Tag = "NON-CANCER/PASSENGER"
        
            if float(Score) < -0.50:
                Tag = "CANCER" 
    else:
        Tag = "TOLERATED"
        
        if float(Score) < -3.00:
            Tag = "DAMAGING"
    
    Predictions.write(
        "\t".join([ str(idx), dbSNP, protein, substitution, Tag, Score, "|".join([ x['description'] for x in Phenotypes ]), "" ]) + "\n"
    ); return True
            
#
if __name__ == '__main__':
    
    #
    # PARSE PROGRAM ARGUMENTS
    #
    
    parser = OptionParser()
    parser.add_option(
                      "-i",
                      dest    = "input",
                      help    = "process dbSNP rs IDs/protein missense mutations from <INPUT>", 
                      metavar = "<INPUT>",
                      default = None
                      )
    parser.add_option(
                      "-o",
                      dest    = "output",
                      help    = "write predictions/phenotype-associations to <OUTPUT>", 
                      metavar = "<OUTPUT>",
                      default = None
                      )
    parser.add_option(
                      "-w",
                      dest    = "weights",
                      help    = "return weighted predictions using pathogenicity weights <WEIGHTS>", 
                      metavar = "<WEIGHTS>",
                      default = "Inherited"
                      )
    parser.add_option(
                      "-p",
                      dest    = "phenotypes",
                      help    = "append domain-phenotype associations for mutations using phenotype ontology <PHENO>", 
                      metavar = "<PHENO>",
                      default = "DO"
                      )
    parser.add_option(
                      "--HTML",
                      dest = "HTML",
                      help = SUPPRESS_HELP,
                      action = "store_true",
                      default = False
                      ) # web-based parameter used to write HTML predictions - hidden from program user(s)

    (options, args) = parser.parse_args()
    
    #
    # AUTHENTICATE PROGRAM PARAMETERS
    #
    
    if not options.input:
        parser.error("No Input File Given (-i parameter)")
    
    if not options.output:
        parser.error("No Output File Given (-o parameter)")
        
    if not options.weights.upper() in [ "UNWEIGHTED", "INHERITED", "CANCER" ]:
        parser.error("Invalid Weighting Scheme")
    
    #
    # INITIALIZE DATABASE CONNECTION/CURSOR
    #
    
    Config = ConfigParser.ConfigParser()
    Config.read("./config.ini")
        
    dbCursor = \
            MySQLdb.connect(
                host     = str(Config.get("DATABASE", "HOST")),
                port     = int(Config.get("DATABASE", "PORT")),
                user     = str(Config.get("DATABASE", "USER")),
                passwd   = str(Config.get("DATABASE", "PASSWD")),
                db       = str(Config.get("DATABASE", "DB")),
                compress = 1
            ).cursor(MySQLdb.cursors.DictCursor)
    
    #
    # PROCESS MUTATION(S)
    #
    
    Predictions = open("/tmp/" + os.path.basename(os.path.splitext(options.input)[0]) + ".tmp", "w")
    Predictions.write("\t".join([ "#", 
                                  "dbSNP ID",
                                  "Protein ID",
                                  "Substitution",
                                  "Prediction",
                                  "Score",
                                  "Domain-Phenotype Association",
                                  "Warning"
                                ]) + "\n")
    
    idx = 0
    for record in open(options.input, "r"):
        if record and not record.startswith("#"):
            record    = record.strip()
            
            try:
                # dbSNP record ...
                if re.compile("^rs\d+$", re.IGNORECASE).match(record):
                    # ... fetch protein consequence(s)
                    dbCursor.execute("select distinct * from Variants where id='" + record + "'")
                    dbRecords = dbCursor.fetchall()
                    
                    if not dbRecords:
                        # no dbSNP mapping
                        idx += 1
                        
                        Predictions.write(
                            "\t".join([ str(idx), record, "", "", "", "", "", "No dbSNP Mapping(s) For Record" ]) + "\n"
                        ); continue
                    
                    # process dbSNP/protein consequence(s) ...
                    for x in dbRecords:
                        idx += 1
                        
                        process_record(x['id'], x['protein'], x['substitution'])
                else:
                    # parse protein/substitution(s) ...
                    Protein       = record.upper().split()[0]
                    Substitutions = record.upper().split()[1].split(",")
                        
                    for x in Substitutions:
                        idx += 1
                        
                        process_record("-", Protein, x)
            #
            except Exception, e:
                idx += 1
                
                Predictions.write(
                    "\t".join([ str(idx), "", "", "", "", "", "", "An Error Occured While Parsing The Record: " + record ]) + "\n"
                )
            #
        #
    #
    
    Predictions.close()
    
    # move predictions to the requested location
    os.system("mv /tmp/" + os.path.basename(os.path.splitext(options.input)[0]) + ".tmp " + options.output)
    
    # 
    # WRITE WEB-BASED PREDICTIONS (IF REQUESTED)
    #
    
    if options.HTML:
        HTM = open(os.path.splitext(options.output)[0] + ".htm", "w")
        HTM.write(
"""
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
        )
        
        for record in open(options.output, "r"):
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
                        
                HTM.write(
"""
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
        )
        
        HTM.write(
"""
    </tbody>
</table>

<script>
    document.getElementById("info").setAttribute("class", "btn btn-primary btn-large pull-right");
    document.getElementById("info").innerHTML = "Download Predictions &raquo;";
    
    clearInterval(Refresh);
</script>
"""
        )
