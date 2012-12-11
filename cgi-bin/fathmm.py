#!/usr/bin/python -u

import os
import re
import math
import ConfigParser

import MySQLdb
import MySQLdb.cursors

from optparse import OptionParser, SUPPRESS_HELP, OptionGroup

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
    # fetch pre-computed sequence record
    #
    
    dbCursor.execute("select a.* from Clusters a, Protein b where a.id=b.cluster and b.name='" + protein + "'")
    SeqRecord = dbCursor.fetchone()
    
    if not SeqRecord:
        # no pre-computed sequence record ...
        if HTM:
            HTM.write(
"""
                    <tr>
                        <td>""" + str(idx) + """</td>
                        <td>""" + str(dbSNP) + """</td>
                        <td>""" + str(protein) + """</td>
                        <td>""" + str(substitution) + """</td>
                        <td colspan="2">&nbsp;</td>
                        <td>No Sequence Record Found</td>
                    </tr>
"""
            )
        TAB.write(
            "\t".join([ str(idx), dbSNP, protein, substitution, "", "", "", "No Sequence Record Found" ]) + "\n"
        ); return False
    
    #
    # authenticate protein/substitution
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
        if HTM:
            HTM.write(
"""
                    <tr>
                        <td>""" + str(idx) + """</td>
                        <td>""" + str(dbSNP) + """</td>
                        <td>""" + str(protein) + """</td>
                        <td>""" + str(substitution) + """</td>
                        <td colspan="2">&nbsp;</td>
                        <td>""" + _ + """</td>
                    </tr>
"""
                )
            TAB.write(
                "\t".join([ str(idx), dbSNP, protein, substitution, "", "", "", _ ]) + "\n"
            ); return False
    
    #
    # create a facade of possible prediction(s)
    #
    
    facade    = []
    
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
                facade.append(dbRecord)
    
    # ... append JackHMMER probabilities
    dbCursor.execute("select a.*, b.accession, b.description from Probabilities a, Library b where a.id=b.id and a.id='" + SeqRecord['id'] + "' and a.position='" + substitution[1:-1] + "'")
    ConRecord = dbCursor.fetchone()
    if ConRecord:
        facade.append(ConRecord)
    
    if not facade:
        # no domain assignments/JackHMMER conservation
        if HTM:
            HTM.write(
"""
                    <tr>
                        <td>""" + str(idx) + """</td>
                        <td>""" + str(dbSNP) + """</td>
                        <td>""" + str(protein) + """</td>
                        <td>""" + str(substitution) + """</td>
                        <td colspan="2">&nbsp;</td>
                        <td>No Prediction Available</td>
                    </tr>
"""
            )
        TAB.write(
            "\t".join([ str(idx), dbSNP, protein, substitution, "", "", "", "No Prediction Available" ]) + "\n"
        ); return False
    
    #
    # derive a prediction using the most informative HMM within the facade
    #
    
    Score      = None
    Phenotypes = None
    
    for x in sorted(facade, key=lambda x:x['information'], reverse=True):
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
    # write our prediction
    #
    
    if not Score:
        # no pathogenicity weights/prediction score
        if HTM:
            HTM.write(
"""
                    <tr>
                        <td>""" + str(idx) + """</td>
                        <td>""" + str(dbSNP) + """</td>
                        <td>""" + str(protein) + """</td>
                        <td>""" + str(substitution) + """</td>
                        <td colspan="2">&nbsp;</td>
                        <td>No Prediction Available</td>
                    </tr>
"""
            )
        TAB.write(
            "\t".join([ str(idx), dbSNP, protein, substitution, "", "", "", "No Prediction Available" ]) + "\n"
        ); return False
    
    Tag = None
    Web = None

    if not options.weights.upper() == "UNWEIGHTED":
        if options.weights.upper() == "INHERITED":
            # "Inherited Disease" predictions ...
            Tag = "TOLERATED"
            Web = "<p class='text-success'>TOLERATED</p>"
        
            if float(Score) < -1.50:
                Tag = "DAMAGING"
                Web = "<p class='text-error'>DAMAGING</p>"    
        else:
            # "Cancer-Associated" predictions ...
            Tag = "OTHER/PASSENGER"
            Web = "<p class='text-success'>OTHER/PASSENGER</p>"
        
            if float(Score) < -0.50:
                Tag = "CANCER"
                Web = "<p class='text-error'>CANCER</p>"    
    else:
        Tag = "TOLERATED"
        Web = "<p class='text-success'>TOLERATED</p>"
        
        if float(Score) < -3.00:
            Tag = "DAMAGING"
            Web = "<p class='text-error'>DAMAGING</p>"
        
    if HTM:
        HTM.write(
"""
                    <tr>
                        <td>""" + str(idx) + """</td>
                        <td>""" + str(dbSNP) + """</td>
                        <td>""" + str(protein) + """</td>
                        <td>""" + str(substitution) + """</td>
                        <td>""" + str(Web) + """</td>
                        <td>""" + str(Score) + """</td>
                        <td><ul>""" + "".join([ "<li>" + x['description'] + "</li>" for x in Phenotypes ]) + """</ul></td>
                    </tr>
"""
        )
    TAB.write(
        "\t".join([ str(idx), dbSNP, protein, substitution, Tag, Score, "|".join([ x['description'] for x in Phenotypes ]), "" ]) + "\n"
    ); return True
            
#
if __name__ == '__main__':
    
    #
    # parse program arguments
    #
    
    parser = OptionParser()
    parser.add_option(
                      "-i",
                      dest    = "batch",
                      help    = "process dbSNP rs IDs/protein missense mutations from <FILE>", 
                      metavar = "<FILE>",
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
                      default = "None"
                      )
    parser.add_option(
                      "-H",
                      dest    = "HTM",
                      help    = SUPPRESS_HELP,
                      action  = "store_true",
                      default = False
                      ) # web-based parameter - hidden from program user(s)
    
    
    group = OptionGroup(parser, "Dangerous Options", "Caution: use these options \n\nat your own risk.  It is believed that some of them bite.")
    
    parser.add_option_group(group)
    
    
    (options, args) = parser.parse_args()
    
    #
    # authenticate program parameters
    #
    
    if not options.batch:
        parser.error("No Mutation Data Given")
        
    if not options.weights.upper() in [ "UNWEIGHTED", "INHERITED", "CANCER" ]:
        parser.error("Invalid Weighting Scheme")
    
    #
    # initialize database connection/cursor
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
    # process mutation(s)
    #
    
    TAB = open(os.path.splitext(options.batch)[0] + ".tab", "w")
    TAB.write("\t".join([ "#", "dbSNP ID", "Protein ID", "Substitution", "Prediction", "Score", "Domain-Phenotype Association", "Warning" ]) + "\n")
    
    HTM = None
    if options.HTM:
        HTM = open(os.path.splitext(options.batch)[0] + ".htm", "w")
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
    
    idx = 0
    for record in open(options.batch, "r"):
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
                        
                        if HTM:
                            HTM.write(
"""
                    <tr>
                        <td>""" + str(idx) + """</td>
                        <td colspan="5">""" + record + """</td>
                        <td>No dbSNP Mapping(s) For Record</td>
                    </tr>
"""
                            )
                        TAB.write(
                            "\t".join([ str(idx), record, "", "", "", "", "", "No dbSNP Mapping(s) For Record" ]) + "\n"
                        ); continue
                    
                    # process dbSNP/protein consequence(s) ...
                    for x in dbRecords:
                        idx += 1; process_record(x['id'], x['protein'], x['substitution'])
                else:
                    # parse protein/substitution(s) ...
                    Protein       = record.upper().split()[0]
                    Substitutions = record.upper().split()[1].split(",")
                        
                    for x in Substitutions:
                        idx += 1; process_record("-", Protein, x)
            #
            except Exception, e:
                # report the exception
                idx += 1
                
                TAB.write(
                    "\t".join([ str(idx), "", "", "", "", "", "", "An Error Occured While Parsing The Record: " + record ]) + "\n"
                )
                if HTM:
                    HTM.write(
"""
                    <tr>
                        <td>""" + str(idx) + """</td>
                        <td colspan="5">&nbsp;</td>
                        <td>An Error Occured While Parsing: """ + record + """</td>
                    </tr>
"""
                    )
            #
        #
    #
    
    if HTM:
        # finish web-based function(s) ...
        HTM.write(
"""
                </tbody>
            </table>
            
            <script>
                document.getElementById("info").setAttribute("class", "btn btn-primary btn-large pull-right");
                document.getElementById("info").innerHTML = "Download Predictions &raquo;";
                document.getElementById("info").href = "../tmp/""" + os.path.basename(os.path.splitext(options.batch)[0]) + """.tab";
                clearInterval(Refresh);
            </script>
"""

        )
