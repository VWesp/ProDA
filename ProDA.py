#--- ProDA: script to locate genes in genomic sequences based on homologous protein sequences using Exonerate and Spaln ---
#loading required modules
import argparse
import yaml
import os
import subprocess
from Bio import SeqIO
import csv
import numpy as np
import traceback
import multiprocessing as mp
from functools import partial
import math
import pandas as pd
import shutil
import re

#option for pandas to display the entire text of a cell in the resulting GFF files
pd.set_option('display.max_colwidth', -1)
#manager and associated variables for multiprocessing
manager = mp.Manager()
#variable to differentiate between the same protein when running the multiprocessing methods
index = manager.Value("i", -1)
#variable for locking variables if multiple instances access them during the multiprocessing methods
lock = manager.Lock()
#variable for displaying the progress of various processes
progress = manager.Value("i", 0)
#variable containing the Fasta file of the current subject as an index
subject_fasta = None

#method for reading the config file
##config_path: path to config file
##return: variable containing the data of the config file as a dictionary
def readConfigFile(config_path):
    config_data = None
    with open(config_path, "r") as config_reader:
        config_data = yaml.safe_load(config_reader)

    return config_data

#method for building a BLAST database of the current subject
##config: dictionary of the config data
##subject: name of the subject in the config file
def buildBLASTDatabase(config, subject):
    #create the folders for the database and the corresponding logs
    db_path = "ProDA/blast_dbs/" + subject + "/"
    log_path = "ProDA/log/blast_dbs/" + subject + "/"
    if(not os.path.exists(db_path)):
        os.makedirs(db_path)

    if(not os.path.exists(log_path)):
        os.makedirs(log_path)

    #access the files in the database folder and check if the 'nhr', 'nin' and 'nsq' files already exist
    #if they do, the database won't be build again
    files_process = subprocess.Popen("ls -p " + db_path + " | grep -v /$",
                    stdout=subprocess.PIPE, shell=True)
    db_files = list(filter(None, files_process.communicate()[0].decode("utf-8").split("\n")))
    if(not (subject+".nhr" in db_files and subject+".nin" in db_files and subject+".nsq" in db_files)):
        print("\t\tBuilding database for subject")
        #bash command for building the database
        db_output = subprocess.Popen("(makeblastdb -in " + config["subjects"][subject] + ""
                    " -dbtype nucl -out " + db_path + subject + ") 2> " + log_path + "makeblastdb.log",
                    stdout=subprocess.PIPE, shell=True)
        res,err = db_output.communicate()
        if(err != None):
            print("\t\t\033[91mError while building database for subject see log file '" + log_path + "'\033[0m")
    else:
        print("\t\t\033[32mFound database files for subject, skipping\033[0m")

#method for running the BLAST command for the current subject and query
##config: dictionary of the config data
##subject: name of the subject in the config file
##query: name of the query in the config file
def runBLASTCommand(config, subject, query):
    #create the folders for the BLAST results and the corresponding logs
    blast_path = "ProDA/blast_results/" + subject + "/"
    log_path = "ProDA/log/blast_results/" + subject + "/"
    if(not os.path.exists(blast_path)):
        os.makedirs(blast_path)

    if(not os.path.exists(log_path)):
        os.makedirs(log_path)

    db_path = "ProDA/blast_dbs/" + subject + "/" + subject
    #check if the required 'fna' file already exist
    #if it does, the hits won't be extracted again
    if(not os.path.exists(blast_path + query + ".fna")):
        print("\t\t\tMatching sequences")
        #check if the resulting BLAST 'hit' file already exist
        #if it does, the BLAST command won't run again
        if(not os.path.exists(blast_path + query + ".hit")):
            print("\t\t\t\tBLAST: evalue=" + str(config["evalue"]))
            #bash command for the running BLAST
            blast_output = subprocess.Popen("(tblastn -query " + config["queries"][query] + ""
                           " -db "+ db_path + " -outfmt '6 qseqid sseqid sstart send evalue'"
                           " -evalue " + str(config["evalue"]) + ""
                           " -num_threads " + str(config["threads"]) + ""
                           " > " + blast_path + query + ".hit" + ")"
                           " 2> " + log_path + query + ".log",
                           stdout=subprocess.PIPE, shell=True)
            res,err = blast_output.communicate()
            if(err != None):
                print("\t\t\t\t\033[91mError for BLAST: subject: " + subject + "\tquery: " + query + ","
                      " see log file '" + log_path + "'\033[0m")
        else:
            print("\t\t\t\t\033[32mFound BLAST file, skipping\033[0m")

        #extract the sequences of the resulting BLAST file
        matchBLASTPositions(blast_path + query + ".hit", config["subjects"][subject],
                            config["left_addendum"], config["right_addendum"], log_path)
        print()
    else:
        print("\t\t\t\t\033[32mFound matching sequence file, skipping\033[0m")

#method for extracting the sequences of each subject for each query ranging from the most left to the most right positions
##hit_path: path to the BLAST file
##subject: name of the subject in the config file
##lad: parameter containing the number nucleotides that should be added to the left side of the extracted sequence
##rad: parameter containing the number nucleotides that should be added to the right side of the extracted sequence
#log: path to the log file
def matchBLASTPositions(hit_path, subject, lad, rad, log):
    print("\t\t\t\tSequence matcher: left=" + str(lad) + "\tright=" + str(rad))
    progress.value = 0
    #get the number of rows in the BLAST file
    length = int(subprocess.Popen("wc -l " + hit_path + " | awk '{ print $1 }'",
                 stdout=subprocess.PIPE, shell=True).communicate()[0].decode("utf-8"))
    print("\t\t\t\t\t" + str(progress.value) + "/" + str(length) + " hits matched", end="\r")
    #check if the BLAST retuned hits
    #if not, skip this process
    if(os.stat(hit_path).st_size != 0):
        log_path = log + hit_path.split("/")[-1].split(".hit")[0] + ".log"
        try:
            #list containing the extracted sequences in fasta format
            fasta_results = []
            #access the BLAST file as a CSV file
            with open(hit_path, "r") as hit_reader:
                hit_content = csv.reader(hit_reader, delimiter="\t")
                current_query = None
                current_subject = None
                current_sequence = None
                #list containing the positive start and end positions of a subject and query
                positions = []
                #list containing the negative (complement strand) start and end positions of a subject and query
                reverse_positions = []
                #loop over the rows of the BLASt file
                for row in hit_content:
                    #if either the query or the subject changes, extract the sequence
                    if(row[0] != current_query or row[1] != current_subject):
                        #extract the positive sequence
                        if(len(positions)):
                            #add the extracted sequence to the result list
                            fasta_results.append(getMatchedSequence(positions, current_sequence, current_query, current_subject, lad, rad))
                            positions = []
                        #extract the negative sequence
                        if(len(reverse_positions)):
                            #add the extracted sequence to the result list
                            fasta_results.append(getMatchedSequence(reverse_positions, current_sequence, current_query, current_subject, lad, rad))
                            reverse_positions = []

                        #get the current query and subject and the sequence of the subject
                        current_query = row[0]
                        current_subject = row[1]
                        current_sequence = str(subject_fasta[current_subject].seq)

                    #get the start and end positions of the local hit
                    start = int(row[2])
                    end = int(row[3])
                    #if the start positions is greater than the end positions (complement strand), add them to the negative list
                    #else, add them to the positive list
                    if(start > end):
                        reverse_positions.append(start)
                        reverse_positions.append(end)
                    else:
                        positions.append(start)
                        positions.append(end)

                    #increase the current by 1 and print it
                    progress.value += 1
                    print("\t\t\t\t\t" + str(progress.value) + "/" + str(length) + " hits matched", end="\r")

                #extract the positive sequence of the last subject and query
                if(len(positions)):
                    #add the extracted sequence to the result list
                    fasta_results.append(getMatchedSequence(positions, current_sequence, current_query, current_subject, lad, rad))
                    positions = []
                #extract the negative sequence of the last subject and query
                if(len(reverse_positions)):
                    #add the extracted sequence to the result list
                    fasta_results.append(getMatchedSequence(reverse_positions, current_sequence, current_query, current_subject, lad, rad))
                    reverse_positions = []

            #write the extracted sequences to the 'fna' file
            with open(hit_path.replace(".hit", ".fna"), "w") as fasta_writer:
                fasta_writer.write("\n".join(fasta_results))

            fasta_results = []
        #catch and print the error if one is thrown
        except:
            print("\t\t\t\t\033[91mError while matching sequences, see log file '" + log_path + "'\033[0m")
            print(str(traceback.format_exc()))
            with open(log_path, "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))
    else:
        with open(hit_path.replace(".hit", ".fna"), "w") as fasta_writer:
            fasta_writer.write("")

#method for extracting the sequence based of the BLAST results
##positions: list containing the local hit positions for a query and subject
##seq: sequence of the current subject
##query: name of the current query
##subject: name of the current subject
##lad: parameter containing the number nucleotides that should be added to the left side of the extracted sequence
##rad: parameter containing the number nucleotides that should be added to the right side of the extracted sequence
##return: extracted sequence in Fasta format
def getMatchedSequence(positions, seq, query, subject, lad, rad):
    #get the most left position of all local hits
    start = min(positions)
    #get the most right position of all local hits
    end = max(positions)
    #add the left addendum to the start position
    #if it is lower than 0, set the start position to 0
    start -= lad
    if(start < 0):
        start = 0
    #add the right addendum to the end position
    #if it is greater than the length of the sequence, set the end position to the length of the sequence
    end += rad
    if(end >= len(seq)):
        end = len(seq) - 1

    header = ">" + query + "::" + subject + "::" + str(start) + "-" + str(end)
    return header + "\n" + seq[start:end]

#method for the running the Exonerate command
##config: dictionary of the config data
##subject: name of the subject in the config file
##query: name of the query in the config file
def runExonerateCommand(config, subject, query):
    #create the folders for the Exonerate results and the corresponding logs
    ex_output = "ProDA/alignment/exonerate/" + subject + "/"
    log_path = "ProDA/log/alignment/exonerate/" + subject + "/"
    if(not os.path.exists(ex_output)):
        os.makedirs(ex_output)

    if(not os.path.exists(log_path)):
        os.makedirs(log_path)

    print("\t\t\tRunning Exonerate: percentage=" + str(config["exonerate_percentage"]) + "\tblosum=" + str(config["blosum"]) + "\trefine=" + config["refine"])
    blast_seq = "ProDA/blast_results/" + subject + "/" + query + ".fna"
    #check if the required 'gff' and 'cds' files already exist
    #if they do, the Exonerate result won't be parsed again
    if(not (os.path.exists(ex_output + query + ".gff") and os.path.exists(ex_output + query + ".cds"))):
        #check if the required 'ryo' file already exists
        #if it does, Exonerate won't be run again
        if(not os.path.exists(ex_output + query + ".ryo")):
            #check if the BLAST file exist
            #if it doesn't, run the BLAST method first
            if(os.path.exists(blast_seq)):
                #check if BLAST returned hits
                #else, skip Exonerate
                if(os.stat(blast_seq).st_size != 0):
                    index.value = -1
                    #access the extracted sequences of the BLAST method
                    matches = list(SeqIO.parse(blast_seq, "fasta"))
                    #set the number of avaiable pools to the number specified in the config file
                    pool = mp.Pool(processes=config["threads"])
                    progress.value = 0
                    #write the file containing the results of the Exonerate command
                    with open(ex_output + query + ".ryo", "w") as ryo_writer:
                        ryo_writer.write("")
                    print("\t\t\t\t" + str(progress.value) + "/" + str(len(matches)) + " sequences processed", end="\r")
                    #start the multiprocessing command for Exonerate for each sequence in the query file and each extracted sequence of the BLAST method
                    pool_map = partial(runExonerateMultiprocessing, query_name=query, queries=config["queries"][query],
                                       percent=config["exonerate_percentage"], blosum=config["blosum"],
                                       refine=config["refine"], output=ex_output, log=log_path, length=len(matches))
                    pool.map_async(pool_map, matches)
                    pool.close()
                    pool.join()
                    print()
                    #parse the Exonerate results
                    print("\t\t\t\tProcessing Exonerate results")
                    parseExonerateResults(ex_output + query + ".ryo")
                else:
                    with open(ex_output + query + ".ryo", "w") as ryo_writer:
                        ryo_writer.write("")
                    with open(ex_output + query + ".gff", "w") as gff_writer:
                        gff_writer.write("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")
                    with open(ex_output + query + ".cds", "w") as cds_writer:
                        cds_writer.write("")
            else:
                print("\t\t\t\t\033[33mRequired files for Exonerate not found, starting Matcher\033[0m")
                #run the BLAST method first and then the Exonerate method
                try:
                    runBLASTCommand(config, subject, query)
                    runExonerateCommand(config, subject, query)
                #catch and print the error if one is thrown
                except:
                    raise Exception(str(traceback.format_exc()))
        else:
            print("\t\t\t\t\033[32mExonerate file found, skipping\033[0m")
            #parse the Exonerate results
            print("\t\t\t\tProcessing Exonerate results")
            parseExonerateResults(ex_output + query + ".ryo")
    else:
        print("\t\t\t\t\033[32mProcessed Exonerate files found, skipping\033[0m")

#method for running the Exonerate parallel for each sequence in the query file and each extracted sequence of the BLAST method
##match: sequence of one extracted BLAST hit in Fasta format
##query_name: name of the query in the config file
##queries: query sequences in Fasta format
##percent: percentage in the config file for the Exonerate command '--percent'
##blosum: number of the BLOSUM matrix in the config file for the Exonerate command '--blosum'
##refine: set options for the Exonerate command '--refine'
##output: path to output directory of the Exonerate results
##log: path to log files of the Exonerate process
##length: number of extracted sequences
def runExonerateMultiprocessing(match, query_name, queries, percent, blosum, refine, output, log, length):
    #set the index for each process
    local_index = None
    with lock:
        index.value += 1
        local_index = index.value

    #get the sequences of the query in the config file as an index
    query_fasta = SeqIO.index(queries, "fasta")
    #get the name of the current query
    query = match.id.split("::")[0]
    #write the extracted sequence to a file
    temp_target = output + query + "_" + str(local_index) + "_target.fna"
    with open(temp_target, "w") as target_writer:
        target_writer.write(">" + match.id + "\n" + str(match.seq))

    #write the current query to a file
    temp_query = output + query + "_" + str(local_index) + "_query.fna"
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query_fasta[query].id + "\n" + str(query_fasta[query].seq))

    temp_output_ryo = output + query + "_" + str(local_index) + ".ryo"
    #run the Exonerate command for the extracted sequence and current correspnding query with the specified commands
    ex_output = subprocess.Popen("(exonerate --model protein2genome --targettype dna --querytype protein"
                " --showtargetgff --showalignment no --showvulgar no --ryo '>%ti\n%tcs'"
                " --refine " + refine + " --proteinsubmat blosum/" + str(config["blosum"]) + ".txt"
                " --percent " + str(config["exonerate_percentage"]) + " --bestn 1"
                " --query " + temp_query + " --target " + temp_target + ""
                " > " + temp_output_ryo + ")"
                " 2> " + log + query + ".log",
                stdout=subprocess.PIPE, shell=True)
    res,err = ex_output.communicate()
    if(err != None):
        print()
        print("\t\t\t\t\033[91mError while running Exonerate, see log file '" + log + query + ".log'\033[0m")

    temp_output_ryo2 = output + query + "_" + str(local_index) + ".ryo2"
    #remove the first line and the last for lines of the local resulting Exonerate file
    os.system("(tail -n +4 " + temp_output_ryo + " | head -n -1) > " + temp_output_ryo2)
    os.remove(temp_target)
    os.remove(temp_query)
    os.remove(temp_output_ryo)
    #write the local results to the overall Exonerate results file and increase the progress by 1 and print it
    with lock:
        with open(temp_output_ryo2, "r") as output_reader:
            with open(output + query_name + ".ryo", "a") as ryo_writer:
                ryo_writer.write(output_reader.read() + "\n")

        os.remove(temp_output_ryo2)
        progress.value += 1

    print("\t\t\t\t" + str(progress.value) + "/" + str(length) + " matches processed", end="\r")

#method for parsing the Exonerate results and converting them to a 'gff' and 'cds' file
##ryo_results: path to Exonerate result file
def parseExonerateResults(ryo_results):
    #check if Exonerate returned results
    #if not, skip the parsing
    if(os.stat(ryo_results).st_size  != 0):
        #list containing the 'gff' format of the Exonerate results
        gff = []
        #list containing the 'cds' format of the Exonerate results
        cds = []
        seq_found = False
        #list containing the exon positions of the current hit
        positions = []
        #access the Exonerate results
        with open(ryo_results, "r") as ryo_reader:
            lines = ryo_reader.readlines()
            #variable containing specific index for each hit
            hit_id = 0
            #loop over each line
            for line in lines:
                if(line.strip()):
                    #if line starts with '#', a new block starts
                    if(line.startswith("#")):
                        seq_found = False
                    #if line starts with '>', the header of the sequence is found
                    elif(line.startswith(">")):
                        header = line.strip().split("::")[0] + "::" + line_splitted[0].split("::")[1] + "::" + str(hit_id)
                        cds.append(header)
                        hit_id += 1
                        seq_found = True
                    #if the line doesn't start with '#' or '>' and the header wasn't found, the subsequent lines are the exons in 'gff' format
                    elif(not seq_found):
                        line_splitted = line.strip().split("\t")
                        query = line_splitted[0].split("::")[0]
                        target = line_splitted[0].split("::")[1]
                        start = int(line_splitted[3]) + int(line_splitted[0].split("::")[2].split("-")[0]) - 1
                        end = int(line_splitted[4]) + int(line_splitted[0].split("::")[2].split("-")[0])
                        #if the current line is the the 'gene' line, check if the positions list contains elements
                        #if it does, a 'gff' block was processed and can be added to the gff list
                        if(line_splitted[2] == "gene"):
                            if(len(positions)):
                                #if the start position of the 'gene' line is greater than the end positions, reverse the positions list to display the 'cds' lines in ascending order
                                if(int(positions[0].split("\t")[5]) > int(positions[-1].split("\t")[5])):
                                    gff.append("\n".join(list(reversed(positions))))
                                else:
                                    gff.append("\n".join(positions))
                            positions = []
                            #add the the 'gff' header and the 'gene' line to the gff list
                            gff.append("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")
                            gff.append(target + "\t" + query + "\texonerate:aln\t" + line_splitted[2] + "\t" + str(start) + "\t" + str(end) + "\t" + line_splitted[6] + "\t" + line_splitted[5] + "\t" + str(hit_id))
                        #if the current line is the the 'cds' line, add it to the positions list
                        elif(line_splitted[2] == "cds"):
                            positions.append(target + "\t" + query + "\texonerate:aln\t" + line_splitted[2] + "\t" + str(start) + "\t" + str(end) + "\t" + line_splitted[6] + "\t" + line_splitted[5] + "\t" + str(hit_id))
                    #if the header was found, the subsequent lines are the sequence
                    elif(seq_found):
                        cds.append(re.sub(r'[^A-Za-z*]', '', line.strip()))
            #if the end of the file is reached, check if the positions list contains elements
            #if it does, the last 'gff' block was processed and can be added to the gff list
            if(len(positions)):
                #if the start position of the 'gene' line is greater than the end positions, reverse the positions list to display the 'cds' lines in ascending order
                if(int(positions[0].split("\t")[5]) > int(positions[-1].split("\t")[5])):
                    gff.append("\n".join(list(reversed(positions))))
                else:
                    gff.append("\n".join(positions))

        #write the gff and cds lists to the corresponding files
        if(not len(gff)):
            gff.append("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")

        with open(ryo_results.replace(".ryo", ".gff"), "w") as gff_writer:
            gff_writer.write("\n".join(gff))

        with open(ryo_results.replace(".ryo", ".cds"), "w") as cds_writer:
            cds_writer.write("\n".join(cds))

        gff = []
        cds = []
        positions = []
    else:
        with open(ryo_results.replace(".ryo", ".gff"), "w") as gff_writer:
            gff_writer.write("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")

        with open(ryo_results.replace(".ryo", ".cds"), "w") as cds_writer:
            cds_writer.write("")

#method for the running the Spaln command
##config: dictionary of the config data
##subject: name of the subject in the config file
##query: name of the query in the config file
def runSpalnCommand(config, subject, query):
    #create the folders for the Spaln results and the corresponding logs
    sp_output = "ProDA/alignment/spaln/" + subject + "/"
    log_path = "ProDA/log/alignment/spaln/" + subject + "/"
    if(not os.path.exists(sp_output)):
        os.makedirs(sp_output)

    if(not os.path.exists(log_path)):
        os.makedirs(log_path)

    print("\t\t\tSpaln: hsp=" + str(config["hsp"]) + "\tstrand=" + str(config["strand"]) + "\tpam=" + str(config["pam"]))
    blast_seq = "ProDA/blast_results/" + subject + "/" + query + ".fna"
    #check if the required 'gff' and 'cds' files already exist
    #if they do, the Spaln result won't be parsed again
    if(not (os.path.exists(sp_output + query + ".gff") and os.path.exists(sp_output + query + ".cds"))):
        #check if the required 'sp' file already exists
        #if it does, Spaln won't be run again
        if(not os.path.exists(sp_output + query + ".sp")):
            #check if the BLAST file exist
            #if it doesn't, run the BLAST method first
            if(os.path.exists(blast_seq)):
                #check if BLAST returned hits
                #else, skip Spaln
                if(os.stat(blast_seq).st_size != 0):
                    index.value = -1
                    #access the extracted sequences of the BLAST method
                    matches = list(SeqIO.parse(blast_seq, "fasta"))
                    #set the number of avaiable pools to the number specified in the config file
                    pool = mp.Pool(processes=config["threads"])
                    progress.value = 0
                     #write the file containing the results of the Spaln command
                    with open(sp_output + query + ".sp", "w") as sp_writer:
                        sp_writer.write("")
                    print("\t\t\t\t" + str(progress.value) + "/" + str(len(matches)) + " sequences processed", end="\r")
                    #start the multiprocessing command for Spaln for each sequence in the query file and each extracted sequence of the BLAST method
                    pool_map = partial(runSpalnMultiprocessing, query_name=query, queries=config["queries"][query], hsp=config["hsp"], strand=config["strand"],
                                       pam=config["pam"], output=sp_output, log=log_path, length=len(matches))
                    pool.map_async(pool_map, matches)
                    pool.close()
                    pool.join()
                    print()
                    #parse the Spaln results
                    print("\t\t\t\tProcessing Spaln results")
                    parseSpalnResults(sp_output + query + ".sp")
                else:
                    with open(sp_output + query + ".sp", "w") as sp_writer:
                        sp_writer.write("")
                    with open(sp_output + query + ".gff", "w") as sp_writer:
                        sp_writer.write("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")
                    with open(sp_output + query + ".cds", "w") as sp_writer:
                        sp_writer.write("")
            else:
                print("\t\t\t\t\033[33mRequired files for Spaln not found, starting Matcher\033[0m")
                #run the BLAST method first and then the Exonerate method
                try:
                    runBLASTCommand(config, subject, query)
                    runSpalnCommand(config, subject, query)
                #catch and print the error if one is thrown
                except:
                    raise Exception(str(traceback.format_exc()))
        else:
            print("\t\t\t\t\033[32mSpaln file found, skipping\033[0m")
            print("\t\t\t\tProcessing Spaln results")
            parseSpalnResults(sp_output + query + ".sp")
    else:
        print("\t\t\t\t\033[32mProcessed Spaln files found, skipping\033[0m")

#method for running the Spaln parallel for each sequence in the query file and each extracted sequence of the BLAST method
##match: sequence of one extracted BLAST hit in Fasta format
##query_name: name of the query in the config file
##queries: query sequences in Fasta format
##hsp: number of the recursive HSP algorithm  in the config file for the Spaln command "-Q"
##strand: number of the specified strand in the config file for the Splan command "-S"
##pam: number of the PAM matrix in the config file for the Spaln commands '-yp' and '-yq'
##output: path to output directory of the Spaln results
##log: path to log files of the Spaln process
##length: number of extracted sequences
def runSpalnMultiprocessing(match, query_name, queries, hsp, strand, pam, output, log, length):
    #set the index for each process
    local_index = None
    with lock:
        index.value += 1
        local_index = index.value

    #get the sequences of the query in the config file as an index
    query_fasta = SeqIO.index(queries, "fasta")
    #get the name of the current query
    query = match.id.split("::")[0]
    #write the extracted sequence to a file
    temp_target = output + query + "_" + str(local_index) + "_target.fna"
    with open(temp_target, "w") as target_writer:
        target_writer.write(">" + match.id + "\n" + str(match.seq))

    #write the current query to a file
    temp_query = output + query + "_" + str(local_index) + "_query.fna"
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query_fasta[query].id + "\n" + str(query_fasta[query].seq))

    temp_output = output + query + "_" + str(local_index) + ".sp"
    #run the Spaln command for the extracted sequence and current correspnding query with the specified commands
    spaln_output = subprocess.Popen("($HOME/spaln2.4.0/bin/spaln -M1 -LS -Q" + str(hsp) + " -S" + str(strand) + " -yp" + str(pam) + " -yq" + str(pam) + ""
                   " -O6 -o" + temp_output + " " + temp_target + " " + temp_query + ")"
                   " 2> " + log + query + ".log",
                   stdout=subprocess.PIPE, shell=True)
    res,err = spaln_output.communicate()
    if(err != None):
        print()
        print("\t\t\t\t\033[91mError while running Spaln, see log file '" + log + query + ".log'\033[0m")

    os.remove(temp_target)
    os.remove(temp_query)
    #write the local results to the overall Spaln results file and increase the progress by 1 and print it
    with lock:
        with open(temp_output, "r") as output_reader:
             with open(output + query_name + ".sp", "a") as sp_writer:
                 sp_writer.write(output_reader.read() + "\n")

        os.remove(temp_output)
        progress.value += 1

    print("\t\t\t\t" + str(progress.value) + "/" + str(length) + " matches processed", end="\r")

#method for parsing the Spaln results and converting them to a 'gff' and 'cds' file
##sp_results: path to Spaln result file
def parseSpalnResults(sp_results):
    #check if Spaln returned results
    #if not, skip the parsing
    if(os.stat(sp_results).st_size != 0):
        #list containing the 'gff' format of the Exonerate results
        gff = []
        #list containing the 'cds' format of the Exonerate results
        cds = []
        #list containing the exon positions of the current hit
        positions = []
        with open(sp_results, "r") as sp_reader:
            lines = sp_reader.readlines()
            target = None
            query = None
            orientation = None
            block_start = None
            #variable containing specific index for each hit
            hit_id = 0
            #loop over each line
            for line in lines:
                #if line starts with ';M', ignore it
                if(line.strip() and not line.startswith(";M")):
                    #if line starts with '>', a new block is reached
                    if(line.startswith(">")):
                        #if the positions list contains elements, add them to the gff list
                        if(len(positions)):
                            for position in positions:
                                    for start_end_pos in position:
                                        #start and end positions are separated by '..'
                                        #Spaln positions are 1-based, so 1 has to be deducted from the start position to be 0-based
                                        start = int(start_end_pos.split("..")[0].strip()) + block_start - 1
                                        end = int(start_end_pos.split("..")[1].strip()) + block_start
                                        gff.append(target + "\t" + query +"\tspaln:aln\tcds\t" + str(start) + "\t" + str(end) + "\t" + orientation + "\t.\t" + str(hit_id))
                            positions = []
                            hit_id += 1

                        #access information of the hit and add them and the 'gff' header to the gff list
                        line_splitted = line.strip().split(" ")
                        query = line_splitted[1].split("::")[0]
                        target = line_splitted[1].split("::")[1]
                        header = ">" + query + "::" + target + "::" + str(hit_id)
                        block_start = int(line_splitted[1].split("::")[2].split("-")[0])
                        start = int(line_splitted[6]) + block_start
                        end = int(line_splitted[8]) + block_start
                        orientation = line_splitted[2]
                        score = line_splitted[-1]
                        gff.append("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")
                        #if the start position is greater then the end position, interchange them
                        if(start > end):
                            temp_start = start
                            start = end
                            end = temp_start
                        start -= 1

                        gff.append(target + "\t" + query + "\tspaln:aln\tgene\t" + str(start) + "\t" + str(end) + "\t" + orientation + "\t" + score + "\t" + str(hit_id))
                        #add the sequence header to the cds list
                        cds.append(header)
                    #if the line starts with ';C', the exon positions are reached
                    elif(line.startswith(";C")):
                        #the first line starts with a 'join' while the subsequent lines do not
                        #exons are separated by ','
                        if("join" in line):
                            positions.append(list(filter(None, line.strip().split("join(")[1].replace(")", "").split(","))))
                        else:
                            positions.append(list(filter(None, line.strip().split(";C ")[1].replace(")", "").split(","))))
                    #all remaining lines contain the sequence
                    else:
                        cds.append(re.sub(r'[^A-Za-z*]', '', line.strip()))

            #if the end of the file is reached, check if the positions list contains elements
            #if it does, the last block was processed and can be added to the gff list
            if(len(positions)):
                for position in positions:
                        for start_end_pos in position:
                            #start and end positions are separated by '..'
                            #Spaln positions are 1-based, so 1 has to be deducted from the start position to be 0-based
                            start = int(start_end_pos.split("..")[0].strip()) + block_start - 1
                            end = int(start_end_pos.split("..")[1].strip()) + block_start
                            gff.append(target + "\t" + query + "\tspaln:aln\tcds\t" + str(start) + "\t" + str(end) + "\t" + orientation + "\t.\t" + str(hit_id))
                positions = []

        #write the gff and cds lists to the corresponding files
        if(not len(gff)):
            gff.append("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")

        with open(sp_results.replace(".sp", ".gff"), "w") as gff_writer:
            gff_writer.write("\n".join(gff))

        with open(sp_results.replace(".sp", ".cds"), "w") as cds_writer:
            cds_writer.write("\n".join(cds))

        gff = []
        cds = []
    else:
        with open(sp_results.replace(".sp", ".gff"), "w") as gff_writer:
            gff_writer.write("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")

        with open(sp_results.replace(".sp", ".cds"), "w") as cds_writer:
            cds_writer.write("")

#method for evaluting the hits generated by Exonerate and Spaln
##config: dictionary of the config data
##algorithm: specify whether Exonerate or Spaln results are evaluated
##subject: name of the subject in the config file
##query: name of the query in the config file
def evaluateAlignment(config, algorithm, subject, query):
    #create the folders for the evaluation results and the corresponding logs
    output = "ProDA/evaluation/" + algorithm + "/" + subject + "/"
    log_path = "ProDA/log/evaluation/" + algorithm + "/" + subject + "/"
    if(not os.path.exists(output)):
        os.makedirs(output)

    if(not os.path.exists(log_path)):
        os.makedirs(log_path)

    print("\t\t\t\tEvaluating results of " + algorithm + ": hssp=" + str(config["hssp_distance"]) + "\tstop-at-stop=" + str(config["stop_at_stop"]) + "\tm-filter=" + str(config["filter_met"]))
    hits_path = "ProDA/alignment/" + algorithm + "/" + subject + "/" + query + ".cds"
    gff = "ProDA/alignment/" + algorithm + "/" + subject + "/" + query + ".gff"
    #check if the required 'stretcher' file exists
    #if it does, the evaluation process won't be run again
    if(not os.path.exists(output + query + ".stretcher")):
        #check if the 'cds' and 'gff' file exist
        #if they don't, run the corresponding algorithm process first
        if(os.path.exists(hits_path) and os.path.exists(gff)):
            #check if the corresponding algorithm returned hits
            #else, skip the evaluation
            if(os.stat(hits_path).st_size != 0):
                index.value = -1
                #access the hit sequences returned by the correspnding algorithm
                hits = list(SeqIO.parse(hits_path, "fasta"))
                #set the number of avaiable pools to the number specified in the config file
                pool = mp.Pool(processes=config["threads"])
                progress.value = 0
                #write the file containing the results of the evaluation
                with open(output + query + ".stretcher", "w") as eval_writer:
                    eval_writer.write("#id\ttarget\tquery\tidentity\tsimilarity\taln_length\thssp_identity\n")
                print("\t\t\t\t\t" + str(progress.value) + "/" + str(len(hits)) + " hits processed", end="\r")
                #start the multiprocessing command for the evaluation for each sequence in the query file and each hit sequence returned by the correspnding algorithm
                pool_map = partial(runStretcherMultiprocessing, query_name=query, queries=config["queries"][query],
                                   distance=config["hssp_distance"],  stop=config["stop_at_stop"], m_fil=config["filter_met"],
                                   output=output, log=log_path, length=len(hits))
                pool.map_async(pool_map, hits)
                pool.close()
                pool.join()
                print()
            else:
                with open(output + query + ".stretcher", "w") as eval_writer:
                    eval_writer.write("#id\ttarget\tquery\tidentity\tsimilarity\taln_length\thssp_identity\n")
        else:
            #run the correspnding algorithm first and then the evaluation
            if(algorithm == "exonerate"):
                print("\t\t\t\t\t\033[33mRequired Exonerate files for evaluation not found, starting Exonerate\033[0m")
                try:
                    runExonerateCommand(config, subject, query)
                #catch and print the error if one is thrown
                except:
                    raise Exception(str(traceback.format_exc()))
            else:
                print("\t\t\t\t\t\033[33mRequired Spaln files for evaluation not found, starting Spaln\033[0m")
                try:
                    runSpalnCommand(config, subject, query)
                #catch and print the error if one is thrown
                except:
                    raise Exception(str(traceback.format_exc()))
            evaluateAlignment(config, algorithm, subject, query)
    else:
        print("\t\t\t\t\t\033[32mEvaluation file found, skipping\033[0m")

#method for running the Stretcher command parallel for each sequence in the query file and each hit sequence returned by Exonerate or Spaln
##hit: sequence of one returned hit
##query_name: name of the query in the config file
##queries: query sequences in Fasta format
##distance: HSSP-distance in the config file
##stop: stop_at_stop number in the config file
##m_fil: filter_met number in the config file
##output: path to output directory of the evaluation results
##log: path to log files of the evaluation process
##length: number of returned hit sequences
def runStretcherMultiprocessing(hit, query_name, queries, distance, stop, m_fil, output, log, length):
    #set the index for each process
    local_index = None
    with lock:
        index.value += 1
        local_index = index.value

    try:
        #get the sequences of the query in the config file as an index
        query_fasta = SeqIO.index(queries, "fasta")
        #get the name of the current query
        query = hit.id.split("::")[0]
        #get the name of the current subject
        subject = hit.id.split("::")[1]
        #get the corresponding index
        id = hit.id.split("::")[2].strip()
        temp_target = output + query + "_" + str(local_index) + "_target.faa"
        translated_sequence = None
        #if the stop paramater is greater than 0, translate the hit sequence for the evaluation only till the first stop codon
        #additionally, trim the sequence
        if(stop > 0):
            translated_sequence = str(trimSequence(hit.seq).translate(to_stop=True))
        else:
            translated_sequence = str(trimSequence(hit.seq).translate())

        #write the current hit sequence to a file
        with open(temp_target, "w") as target_writer:
            target_writer.write(">" + hit.id + "\n" + translated_sequence)

        temp_query = output + query + "_" + str(local_index) + "_query.faa"
        #write the current query sequence to a file
        with open(temp_query, "w") as query_writer:
            query_writer.write(">" + query_fasta[query].id + "\n" + str(query_fasta[query].seq))

        #variable containing the results of the evaluation
        st_result = ""
        #check if the m_filter is greater than 0
        #if it is, the hit sequence has to start with a methionine
        #else, it doesn't matter
        if(m_fil <= 0 or translated_sequence.startswith("M")):
            temp_output = output + query + "_" + str(local_index) + ".stretcher"
            #run the Stretcher command
            os.system("(stretcher -asequence " + temp_query + " -sprotein1"
                      " -bsequence " + temp_target + " -sprotein2 -auto -stdout > " + temp_output + ")"
                      " 2> " + log + query + ".log")
            al_length = None
            identity = None
            similarity = None
            #access the Stretcher result file
            with open(temp_output, "r") as output_reader:
                    content = output_reader.readlines()
                    #loop over each line and save the length of the alignment, the identity percentage and the similarity percentage
                    for line in content:
                        if(line):
                            if(line.startswith("# Length")):
                                al_length = int(line.split(" ")[2].strip())
                            if(line.startswith("# Identity")):
                                identity = float(line.split("(")[1][:-3].strip())
                            if(line.startswith("# Similarity")):
                                similarity = float(line.split("(")[1][:-3].strip())
                                break

            #calculate the HSSP-distance of the resulting alignment
            hssp_identity = calculateHSSPIdentity(al_length, distance)
            #if the identity percentage is greater then the HSSP identity, save the results
            #else, discard them
            if(identity >= hssp_identity):
                st_result = "\t".join([id, subject, query, str(identity), str(similarity), str(al_length), str(hssp_identity)])
            os.remove(temp_output)

        os.remove(temp_target)
        os.remove(temp_query)
        #write the local results to the overall evaluation results file and increase the progress by 1 and print it
        with lock:
            with open(output + query_name + ".stretcher", "a") as eval_writer:
                eval_writer.write(st_result + "\n")

            progress.value += 1

        print("\t\t\t\t\t" + str(progress.value) + "/" + str(length) + " hits processed", end="\r")
    #catch and print the error if one is thrown
    except:
        print("\t\t\t\t\t\033[91mError while evaluating, see log file '" + log + query + ".log'\033[0m")
        print(str(traceback.format_exc()))
        with open(log + query + ".log", "w") as log_writer:
            log_writer.write(str(traceback.format_exc()))

#method for trimming a sequence so 'seq_len % 3 == 0'
##seq: DNA sequence
##return: trimmed sequence
def trimSequence(seq):
    return seq[:len(seq)-len(seq)%3]

#method for calculating the HSSP identity of an alignment based on its length
##lenght: length of the alignment
##distance: distance added to the calculated HSSP identity
##return: HSSP identity
def calculateHSSPIdentity(length, distance):
    #if the alignment length is greater than 418, return 19.5 (plus the distance)
    #else, calculate the HSSP identity (plus the distance)
    #if the resulting identity is greater than 100, return 100
    if(length >= 418):
        if(distance + 19.5 > 100):
            return 100
        else:
            return distance + 19.5
    else:
        hssp_distance = distance + 480 * length**(-0.32 * (1 + math.exp(-length / 1000)))
        if(hssp_distance > 100):
            return 100
        else:
            return hssp_distance

#method for merging the evaluated Exonerate and Spaln results
##config: dictionary of the config data
##subject: name of the subject in the config file
##query: name of the query in the config file
def mergeResults(config, subject, query):
    print("\t\t\tMerging duplicates: overlap percentage=" + str(config["overlap_percentage"]))
    try:
        #create the folders for the final results and the corresponding logs
        output = "ProDA/results/" + subject + "/" + query + "/"
        log_path = "ProDA/log/results/" + subject + "/" + query + "/"
        if(not os.path.exists(output)):
            os.makedirs(output)

        if(not os.path.exists(log_path)):
            os.makedirs(log_path)

        #list containing the merged 'gff' result
        gff_list = ["#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid"]
        #list containing the merged 'cds' result
        cds_list = []
        #path to the 'gff' and 'stretcher' files of Exonerate
        ex_gff = "ProDA/alignment/exonerate/" + subject + "/" + query + ".gff"
        ex_stretcher = "ProDA/evaluation/exonerate/" + subject + "/" + query + ".stretcher"
        #list containing processed IDs of the Exonerate files
        ex_id_list = []
        #path to the 'gff' and 'stretcher' files of Spaln
        sp_gff = "ProDA/alignment/spaln/" + subject + "/" + query + ".gff"
        sp_stretcher = "ProDA/evaluation/spaln/" + subject + "/" + query + ".stretcher"
        #list containing processed IDs of the Spaln files
        sp_id_list = []
        #variable containing specific index for each hit
        hit_id = 0
        #check if the merged 'gff', 'nuc', 'cds' and 'pep' files already exist
        #if they do, skip the merging
        if(not (os.path.exists(output + query + ".gff") or os.path.exists(output + query + ".nuc") or os.path.exists(output + query + ".cds") or os.path.exists(output + query + ".pep"))):
            #check if the required 'gff' and 'stretcher' files of the Exonerate and Spaln process exist
            #if not, run the Exonerate, Spaln and evaluation processes first and then merge
            if(os.path.exists(ex_gff) and os.path.exists(ex_stretcher) and os.path.exists(sp_gff) and os.path.exists(sp_stretcher)):
                #get the 'cds' sequences returned by Spaln
                sp_sequences = SeqIO.index("ProDA/alignment/spaln/" + subject + "/" + query + ".cds", "fasta")
                #get the 'gff' file returned by Spaln
                sp_gff_df = pd.read_csv(sp_gff, sep="\t", header=0).astype(str)
                #get the evaluation file of the Spaln results
                sp_eval_df = pd.read_csv(sp_stretcher, sep="\t", header=0).astype(str)
                #get the 'cds' sequences returned by Exonerate
                ex_sequences = SeqIO.index("ProDA/alignment/exonerate/" + subject + "/" + query + ".cds", "fasta")
                #get the 'gff' file returned by Exonerate
                ex_gff_df = pd.read_csv(ex_gff, sep="\t", header=0).astype(str)
                #get the evaluation file of the Exonerate results
                ex_eval_df = pd.read_csv(ex_stretcher, sep="\t", header=0).astype(str)
                progress.value = 0
                length = len(sp_eval_df) * 2 + len(ex_eval_df)
                print("\t\t\t\t" + str(progress.value) + "/" + str(length) + " hits processed", end="\r")
                #loop over the lines of the Spaln evaluation file
                for sp_index,sp_row in sp_eval_df.iterrows():
                    sp_header = sp_row["query"] + "::" + sp_row["target"] + "::" + sp_row["#id"]
                    #get the Spaln 'cds' sequence of the current query, subject and index
                    sp_seq = str(sp_sequences[sp_header].seq)
                    #get the 'gene' line in the Spaln 'gff' file of the current query, subject and index
                    sp_gene_rows = sp_gff_df.loc[(sp_gff_df["type"] == "gene") & (sp_gff_df["#subject"] == sp_row["target"]) & (sp_gff_df["query"] == sp_row["query"])]
                    sp_id_row = sp_gene_rows.loc[sp_gene_rows["id"] == sp_row["#id"]].iloc[0]
                    #get the all 'gene' line in the Exonerate 'gff' file of the current query and subject
                    ex_gene_rows = ex_gff_df.loc[(ex_gff_df["type"] == "gene") & (ex_gff_df["#subject"] == sp_row["target"]) & (ex_gff_df["query"] == sp_row["query"])]
                    #get the all line in the Exonerate evaluation file of the current query and subject
                    ex_eval_rows = ex_eval_df.loc[(ex_eval_df["target"] == sp_row["target"]) & (ex_eval_df["query"] == sp_row["query"])]
                    #loop over the retrieved lines of the Exonerate evaluation file
                    for ex_index,ex_row in ex_eval_rows.iterrows():
                        if(not ex_row["#id"] in ex_id_list):
                            ex_header = ex_row["query"] + "::" + ex_row["target"] + "::" + ex_row["#id"]
                            #get the Exonerate 'cds' sequence of the current query, subject and index
                            ex_seq = str(ex_sequences[ex_header].seq)
                            #get the 'gene' line in the Exonerate 'gff' file of the index
                            ex_id_row = ex_gene_rows.loc[ex_gene_rows["id"] == ex_row["#id"]].iloc[0]
                            #calculate how much the start and end positions of both hits overlap
                            #if it is higher than the value specified in the config file ('overlap_percentage'), check further
                            #else, continue to the next hit
                            if(calculateOverlapPercentage(int(sp_id_row["start"]), int(sp_id_row["end"]), int(ex_id_row["start"]), int(ex_id_row["end"])) >= config["overlap_percentage"]):
                                #add the IDs of the current Spaln and Exonerate lines to their corresponding lists
                                sp_id_list.append(sp_row["#id"])
                                ex_id_list.append(ex_row["#id"])
                                id_rows = None
                                header = None
                                query_name = None
                                subject_name = None
                                seq = None
                                #check if the Spaln and Exonerate sequences start with an 'ATG' or not
                                #if they do, check further
                                #if only one starts with an 'ATG', the other one will be discarded
                                if((sp_seq.startswith("ATG") and ex_seq.startswith("ATG"))
                                   or not (sp_seq.startswith("ATG") or ex_seq.startswith("ATG"))):
                                   #check what hit has the higher score
                                   #discard the other one
                                   if(float(sp_id_row["score"]) >= float(ex_id_row["score"])):
                                       #get the entire 'gff' block of the hit (specified by the index)
                                       id_rows = sp_gff_df.loc[sp_gff_df["id"] == str(sp_row["#id"])]
                                       header = "##" +  sp_row["target"] + "\t" + sp_row["query"] + "\tid:" + sp_row["identity"] + "%\tsim:" + sp_row["similarity"] + "%\t" + str(hit_id)
                                       query_name = sp_row["query"]
                                       subject_name = sp_row["target"]
                                       seq = sp_seq
                                   else:
                                        #get the entire 'gff' block of the hit (specified by the index)
                                       id_rows = ex_gff_df.loc[ex_gff_df["id"] == ex_row["#id"]]
                                       header = "##" +  ex_row["target"] + "\t" + ex_row["query"] + "\tid:" + ex_row["identity"] + "%\tsim:" + ex_row["similarity"] + "%\t" + str(hit_id)
                                       query_name = ex_row["query"]
                                       subject_name = ex_row["target"]
                                       seq = ex_seq
                                elif(sp_seq.startswith("ATG")):
                                     #get the entire 'gff' block of the hit (specified by the index)
                                    id_rows = sp_gff_df.loc[sp_gff_df["id"] == sp_row["#id"]]
                                    header = "##" +  sp_row["target"] + "\t" + sp_row["query"] + "\tid:" + sp_row["identity"] + "%\tsim:" + sp_row["similarity"] + "%\t" + str(hit_id)
                                    query_name = sp_row["query"]
                                    subject_name = sp_row["target"]
                                    seq = sp_seq
                                else:
                                     #get the entire 'gff' block of the hit (specified by the index)
                                    id_rows = ex_gff_df.loc[ex_gff_df["id"] == ex_row["#id"]]
                                    header = "##" +  ex_row["target"] + "\t" + ex_row["query"] + "\tid:" + ex_row["identity"] + "%\tsim:" + ex_row["similarity"] + "%\t" + str(hit_id)
                                    query_name = ex_row["query"]
                                    subject_name = ex_row["target"]
                                    seq = ex_seq

                                #add the 'gff' block and 'cds' sequence of the winner to the gff and cds lists
                                id_rows["id"] = str(hit_id)
                                #change the accessed 'gff' block to a string without indices and header
                                rows = id_rows.to_string(header=False, index=False).split("\n")
                                #remove all empty columns
                                formatted_rows = ["\t".join(list(filter(None, row.split(" ")))) for row in rows]
                                gff_list.append(header + "\n" + "\n".join(formatted_rows))
                                cds_list.append(">" + query_name + "::" + subject_name + "::" + str(hit_id) + "\n" + sp_seq)
                                hit_id += 1
                                break

                    #increase the progress by 1 and print it
                    progress.value += 1
                    print("\t\t\t\t" + str(progress.value) + "/" + str(length) + " hits processed", end="\r")

                #loop over the lines of the Spaln evaluation file but ignore hits that were already processed (sp_id_list)
                for sp_index,sp_row in sp_eval_df.iterrows():
                    if(not sp_row["#id"] in sp_id_list):
                        #get the entire 'gff' block of the hit (specified by the index)
                        id_rows = sp_gff_df.loc[sp_gff_df["id"] == sp_row["#id"]]
                        id_rows["id"] = str(hit_id)
                        header = "##" +  sp_row["target"] + "\t" + sp_row["query"] + "\tid:" + sp_row["identity"] + "%\tsim:" + sp_row["similarity"] + "%\t" + str(hit_id)
                        sp_header = sp_row["query"] + "::" + sp_row["target"] + "::" + sp_row["#id"]
                        #change the accessed 'gff' block to a string without indices and header
                        rows = id_rows.to_string(header=False, index=False).split("\n")
                        #remove all empty columns
                        formatted_rows = ["\t".join(list(filter(None, row.split(" ")))) for row in rows]
                        gff_list.append(header + "\n" + "\n".join(formatted_rows))
                        cds_list.append(">" + sp_row["query"] + "::" + sp_row["target"] + "::" + str(hit_id) + "\n" + str(sp_sequences[sp_header].seq))
                        hit_id += 1

                    #increase the progress by 1 and print it
                    progress.value += 1
                    print("\t\t\t\t" + str(progress.value) + "/" + str(length) + " hits processed", end="\r")

                #loop over the lines of the Exonerate evaluation file but ignore hits that were already processed (ex_id_list)
                for ex_index,ex_row in ex_eval_df.iterrows():
                    if(not ex_row["#id"] in ex_id_list):
                        #get the entire 'gff' block of the hit (specified by the index)
                        id_rows = ex_gff_df.loc[ex_gff_df["id"] == ex_row["#id"]]
                        id_rows["id"] = str(hit_id)
                        header = "##" +  ex_row["target"] + "\t" + ex_row["query"] + "\tid:" + ex_row["identity"] + "%\tsim:" + ex_row["similarity"] + "%\t" + str(hit_id)
                        ex_header = ex_row["query"] + "::" + ex_row["target"] + "::" + ex_row["#id"]
                        #change the accessed 'gff' block to a string without indices and header
                        rows = id_rows.to_string(header=False, index=False).split("\n")
                        #remove all empty columns
                        formatted_rows = ["\t".join(list(filter(None, row.split(" ")))) for row in rows]
                        gff_list.append(header + "\n" + "\n".join(formatted_rows))
                        cds_list.append(">" + ex_row["query"] + "::" + ex_row["target"] + "::" + str(hit_id) + "\n" + str(ex_sequences[ex_header].seq))
                        hit_id += 1

                    #increase the progress by 1 and print it
                    progress.value += 1
                    print("\t\t\t\t" + str(progress.value) + "/" + str(length) + " hits processed", end="\r")

                #write the gff list to a file
                with open(output + query + ".gff", "w") as gff_writer:
                    gff_writer.write("\n".join(gff_list))

                progress.value = 0
                length = int(subprocess.Popen("wc -l " + output + query + ".gff" + " | awk '{ print $1 }'",
                             stdout=subprocess.PIPE, shell=True).communicate()[0].decode("utf-8"))
                print("\n\t\t\t\t" + str(progress.value) + "/" + str(length) + " GFF lines processed", end="\r")
                nucl = []
                #access the gff file
                with open(output + query + ".gff", "r") as gff_reader:
                    content = csv.reader(gff_reader, delimiter="\t")
                    next(content)
                    target = None
                    protein = None
                    id = None
                    #loop over each row and use the start and end positions of the 'gene' line to get the nucleotide sequence in Fasta format
                    for row in content:
                        if(row[0].startswith("##")):
                            target = row[0][2:]
                            protein = row[1]
                            id = row[4]
                        elif(row[3] == "gene"):
                            start = int(row[4])
                            end = int(row[5])
                            nucl.append(">" + protein + "::" + target + "::" + id + "\n" + str(subject_fasta[target].seq)[start:end])

                        #increase the progress by 1 and print it
                        progress.value += 1
                        print("\t\t\t\t" + str(progress.value) + "/" + str(length) + " GFF lines processed", end="\r")

                #write the nucleotide Fasta to a file
                with open(output + query + ".nuc", "w") as nuc_writer:
                    nuc_writer.write("\n".join(nucl))

                #write the cds Fasta to a file
                with open(output + query + ".cds", "w") as cds_writer:
                    cds_writer.write("\n".join(cds_list))

                #loop over the cds sequences in Fasta format and translate them
                translated_cds = []
                for fasta in SeqIO.parse(output + query + ".cds", "fasta"):
                    translated_cds.append(">" + fasta.id + "\n" + str(fasta.seq.translate()))

                #write the translated cds sequences in Fasta format to a file
                with open(output + query + ".pep", "w") as pep_writer:
                    pep_writer.write("\n".join(translated_cds))

                gff_list = []
                nucl = []
                cds_list = []
                translated_cds = []
                ex_id_list = []
                sp_id_list = []
            else:
                #run the Exonerate process and evaluation first and then merge
                if(not os.path.exists(ex_gff) or not os.path.exists(ex_stretcher)):
                    print("\t\t\t\t\033[33mRequired Exonerate files for Merger not found, starting Exonerate\033[0m")
                    try:
                        runExonerateCommand(config, subject, query)
                        evaluateAlignment(config, "exonerate", subject, query)
                    except:
                        raise Exception(str(traceback.format_exc()))

                #run the Spaln process and evaluation first and then merge
                if(not os.path.exists(sp_gff) or not os.path.exists(sp_stretcher)):
                    print("\t\t\t\t\033[33mRequired Spaln files for Merger not found, starting Spaln\033[0m")
                    try:
                        runSpalnCommand(config, subject, query)
                        evaluateAlignment(config, "spaln", subject, query)
                    except:
                        raise Exception(str(traceback.format_exc()))

                mergeResults(config, subject, query)
        else:
            print("\t\t\t\t\033[32mMerged result files found, skipping\033[0m")
    #catch and print the error if one is thrown
    except:
        print("\t\t\t\t\033[91mError while merging results, see log file '" + log_path + query + ".log'\033[0m")
        print(str(traceback.format_exc()))
        with open(log_path + query + ".log", "w") as log_writer:
            log_writer.write(str(traceback.format_exc()))

#calculate the overlap percentage between two hits based on their start and end positions
##start1: start position of the 1st hit
##end1: end position of the 1st hit
##start2: start position of the 2nd hit
##end2: end position of the 2nd hit
##return: overlap percentage
def calculateOverlapPercentage(start1, end1, start2, end2):
    #check if the two hits overlap at all
    #if they don't, return 0
    #else, calculate the overlap percentage
    if(end1 < start2 or end2 < start1):
        return 0
    else:
        return ((min(end1, end2) - max(start1, start2)) / min(end1-start1, end2-start2)) * 100

#method for extracting the best hit for each query in the merged gff file
##config: dictionary of the config data
##subject: name of the subject in the config file
##query: name of the query in the config file
def extractBestHits(config, subject, query):
    print("\n\t\t\tExtracting best hits: score=" + config["best"] + "\tpre-introns=" + str(config["pre_introns"]))
    output = "ProDA/results/" + subject + "/" + query + "/" + query
    #access the merged gff file
    with open(output + ".gff", "r") as gff_reader:
        content = csv.reader(gff_reader, delimiter="\t")
        next(content)
        header = None
        protein = None
        identity = -1
        similarity = -1
        score = -1
        index = None
        introns = -1
        #dictionary containing the best hit of each query
        best_hit = {}
        #list containing the 'gene' and 'cds' lines of the last block
        current_cds = []
        #loop over each row
        for row in content:
            #if the lien startswith '##' a new block is reached
            if(row[0].startswith("##")):
                #check if the current_cds list contains element
                #if it does, a block was processed and check further
                #else, no block was processed yet
                if(len(current_cds)):
                    #check if the processed query is in the best_hit dictionary
                    #if not, this query is for now the best hit and add it
                    #else check further
                    if(not protein in best_hit):
                        best_hit[protein] = {}
                        best_hit[protein]["i"] = identity
                        best_hit[protein]["s"] = similarity
                        best_hit[protein]["a"] = score
                        best_hit[protein]["c"] =  contig
                        best_hit[protein]["ix"] = index
                        best_hit[protein]["in"] = introns
                        best_hit[protein]["g"] = [header + "\n" + "\n".join(current_cds)]
                    else:
                        #variable checking whether hits with introns are preferred over ones without introns and whether the new query contains introns but not the current best one
                        introns_found = config["pre_introns"] > 0 and introns > 0 and best_hit[protein]["in"] == 0
                        #variable checking whether hits with introns are preferred over ones without introns and whether the the current best one contains introns but not the new query
                        no_introns_found = config["pre_introns"] > 0 and introns == 0 and best_hit[protein]["in"] > 0
                        #variable checking if the new query is better than the current best one based identity, similarity or alignment score
                        score_found = (config["best"] == "i" and identity > best_hit[protein]["i"]) or (config["best"] == "s" and similarity > best_hit[protein]["s"]) or (config["best"] == "a" and score > best_hit[protein]["a"])
                        #check all three conditions
                        #if introns are preferred, check if the new query contains introns but the not the best one or if the new query and the best one both contain introns or not and if the taken score of the new query is better than the best one
                        #if introns are not preferred, only check if the taken score of the new query is better than the best one
                        #else, ignore the new query
                        if(introns_found or (not no_introns_found and score_found)):
                            best_hit[protein] = {}
                            best_hit[protein]["i"] = identity
                            best_hit[protein]["s"] = similarity
                            best_hit[protein]["a"] = score
                            best_hit[protein]["c"] =  contig
                            best_hit[protein]["ix"] = index
                            best_hit[protein]["in"] = introns
                            best_hit[protein]["g"] = [header + "\n" + "\n".join(current_cds)]

                #get the 'gff' header line and the associated information
                header = "\t".join(row)
                contig = row[0][2:]
                protein = row[1]
                identity = float(row[2].split("id:")[1][:-1])
                similarity = float(row[3].split("sim:")[1][:-1])
                introns = -1
                current_cds = []
            else:
                #get the 'gene' and 'cds' lines
                #if it is 'cds' line, increase the number of introns by 1
                current_cds.append("\t".join(row))
                if(row[3] == "gene"):
                    score = float(row[7])
                    index = row[8]
                elif(row[3] == "cds"):
                    introns += 1

        #if the end of the line is reached, check if the current_cds list contains element
        if(len(current_cds)):
            #check if the last processed query is in the best_hit dictionary
            #if not, this query is for now the best hit and add it
            #else check further
            if(not protein in best_hit):
                best_hit[protein] = {}
                best_hit[protein]["i"] = identity
                best_hit[protein]["s"] = similarity
                best_hit[protein]["a"] = score
                best_hit[protein]["c"] =  contig
                best_hit[protein]["ix"] = index
                best_hit[protein]["in"] = introns
                best_hit[protein]["g"] = [header + "\n" + "\n".join(current_cds)]
            else:
                #variable checking whether hits with introns are preferred over ones without introns and whether the new query contains introns but not the current best one
                introns_found = config["pre_introns"] > 0 and introns > 0 and best_hit[protein]["in"] == 0
                #variable checking whether hits with introns are preferred over ones without introns and whether the the current best one contains introns but not the new query
                no_introns_found = config["pre_introns"] > 0 and introns == 0 and best_hit[protein]["in"] > 0
                #variable checking if the new query is better than the current best one based identity, similarity or alignment score
                score_found = (config["best"] == "i" and identity > best_hit[protein]["i"]) or (config["best"] == "s" and similarity > best_hit[protein]["s"]) or (config["best"] == "a" and score > best_hit[protein]["a"])
                #check all three conditions
                #if introns are preferred, check if the new query contains introns but the not the best one or if the new query and the best one both contain introns or not and if the taken score of the new query is better than the best one
                #if introns are not preferred, only check if the taken score of the new query is better than the best one
                #else, ignore the new query
                if(introns_found or (not no_introns_found and score_found)):
                    best_hit[protein] = {}
                    best_hit[protein]["i"] = identity
                    best_hit[protein]["s"] = similarity
                    best_hit[protein]["a"] = score
                    best_hit[protein]["c"] =  contig
                    best_hit[protein]["ix"] = index
                    best_hit[protein]["in"] = introns
                    best_hit[protein]["g"] = [header + "\n" + "\n".join(current_cds)]

            current_cds = []

        #write the the best hits of the gff file to a new gff file
        with open(output + "_best_hits.gff", "w") as gff_writer:
            gff_writer.write("")

        with open(output + "_best_hits.gff", "a") as gff_writer:
            gff_writer.write("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid\n")
            for hit in best_hit:
                gff_writer.write("\n".join(best_hit[hit]["g"]) + "\n")

        #loop over the associated files and extract the sequences (nucleotide, cds and peptide) of the best hits and write them to a new file
        for extension in [".nuc", ".cds", ".pep"]:
            fasta = SeqIO.index(output + extension, "fasta")
            with open(output + "_best_hits" + extension, "w") as fasta_writer:
                fasta_writer.write("")
            for protein in best_hit:
                index = protein + "::" + best_hit[protein]["c"] + "::" + best_hit[protein]["ix"]
                with open(output + "_best_hits" + extension, "a") as fasta_writer:
                    fasta_writer.write(">" + fasta[index].id + "\n" + str(fasta[index].seq) + "\n")

#method for visualizing the best hits of the results
##config: dictionary of the config data
##subject: name of the subject in the config file
##query: name of the query in the config file
def visualizeResults(config, subject, query):
    print("\t\t\tVisualizing results: properties file: " + config["properties"])
    protein = "ProDA/results/" + subject + "/" + query + "/" + query + "_best_hits.pep"
    output = "ProDA/results/" + subject + "/" + query + "/alignment/"
    #check if the visualization folder already exists
    #if it does, skip the visualization
    if(not os.path.exists(output)):
        #check if the required protein file of the best hits exist
        #if it doesn't, run the merging and hit extracting processes first
        if(os.path.exists(protein)):
            #create the folders for the visualizing and the corresponding logs
            os.makedirs(output)
            log_path = "ProDA/log/results/" + subject + "/" + query + "/alignment/"
            if(not os.path.exists(log_path)):
                os.makedirs(log_path)

            #access the protein sequences of the best hits
            protein_fasta = list(SeqIO.parse(protein, "fasta"))
            #set the number of avaiable pools to the number specified in the config file
            pool = mp.Pool(processes=config["threads"])
            progress.value = 0
            print("\t\t\t\t" + str(progress.value) + "/" + str(len(protein_fasta)) + " hits visualized", end="\r")
            #start the multiprocessing command for the visualization for each sequence in the protein file of the best hits
            pool_map = partial(runJalViewMultiprocessing, queries=config["queries"][query],
                               properties=config["properties"], output=output, log=log_path, length=len(protein_fasta))
            pool.map_async(pool_map, protein_fasta)
            pool.close()
            pool.join()
            print()
        else:
            print("\t\t\t\t\033[33mRequired files for Visualizer not found, starting Merger\033[0m")
            #merge the results of Exonerate and Spaln and extract the best hits first and then visualize the results
            try:
                mergeResults(config, subject, query)
                extractBestHits(config, subject, query)
                visualizeResults(config, subject, query)
            #catch and print the error if one is thrown
            except:
                raise Exception(str(traceback.format_exc()))
    else:
        print("\t\t\t\t\033[32mVisualizing folder found, skipping\033[0m")

#method for visualizing the results parallel with Stretcher and JalView
##match: sequence of one extracted BLAST hit in Fasta format
#properties: path to properties file in the config file for JalView
##output: path to output directory of the Exonerate results
##log: path to log files of the Exonerate process
##length: number of extracted sequences
def runJalViewMultiprocessing(hit, queries, properties, output, log, length):
    #get the sequences of the query in the config file as an index
    query_fasta = SeqIO.index(queries, "fasta")
    #get the name of the current query
    query = hit.id.split("::")[0]
    #get the index of the current query
    id = hit.id.split("::")[2]
    #write the best hit protein sequence to a file
    temp_target = output + query + "_target.fna"
    with open(temp_target, "w") as target_writer:
        target_writer.write(">" + hit.id + "\n" + str(hit.seq))

    #write the current query to a file
    temp_query = output + query + "_query.fna"
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query_fasta[query].id + "\n" + str(query_fasta[query].seq))

    clustal = output + query + "_" + id + ".clustal"
    #run the Stretcher command to build an alignment in clustal format
    os.system("(stretcher -asequence " + temp_query + " -sprotein1 -bsequence " + temp_target + ""
              " -sprotein1 -auto -aformat clustal -stdout > " + clustal + ") 2> " + log + query + "_" + id + ".cl_log")
    png = output + query + "_" + id + ".png"
    #run the JalView command to visualize the alignment in PNG format
    jalview_output = subprocess.Popen("(jalview -nodisplay -props " + properties + " -colour clustal"
                     " -open " + clustal + " -png " + png + ") 2> " + log + query + "_" + id + ".jv_log",
                     stdout=subprocess.PIPE, shell=True)
    res,err = jalview_output.communicate()
    if(err != None):
        print("\t\t\t\t\033[91mError while visualizing results, see log file '" + log + query + "_" + id + ".jv_log'\033[0m")

    os.remove(temp_target)
    os.remove(temp_query)
    #increase the progress by 1 and print it
    with lock:
        progress.value += 1

    print("\t\t\t\t" + str(progress.value) + "/" + str(length) + " hits visualized", end="\r")

#command line interface to parse the config file (default: same directory as the script) and force a rerun
parser = argparse.ArgumentParser(description="ProDA (Protein-DNA-Aligner): Search for protein sequences in DNA sequences with the Exonerate and Spaln alignment tools and analyze the results.")
parser.add_argument("-c", "--config", help="Path to config.yaml file", default="config.yaml")
parser.add_argument("-f", "--force", help="Force a rerun of the procedure involving all steps", action="store_true")
args = parser.parse_args()
#check if the config file exists
#if not, throw an error
if(not os.path.exists(args.config)):
    print("\033[91mError: Unable to find/access config file named '" + args.config + "'\033[0m.")
else:
    #check if the 'force' option is set
    #if it is, delete the output folder if it exists
    if(args.force):
        print("Deleting ProDA files")
        if(os.path.exists("ProDA")):
            shutil.rmtree("ProDA")

        if(os.path.exists("ProDA.zip")):
            os.remove("ProDA.zip")

    #create the 'ProDA' output folder
    if(not os.path.exists("ProDA")):
        os.makedirs("ProDA")

    #read the config file as a dictionary
    config = readConfigFile(args.config)
    print("Running ProDA")
    try:
        #loop over each subject in the config file
        for subject in config["subjects"]:
            print("\tSubject: " + subject)
            print("\t\tLoading subject file")
            #replace the symbols ' ' and '|' with a '_' in the headers of the subject file
            os.system("sed -i 's/[ |]/_/g' " + config["subjects"][subject])
            #access the subject sequences as an index
            subject_fasta = SeqIO.index(config["subjects"][subject], "fasta")
            #build a BLASt database of the subject sequences
            buildBLASTDatabase(config, subject)
            #loop over each query in the config file
            for query in config["queries"]:
                print("\t\tQuery: " + query)
                #start the 'ProDA' process and visualize the results
                visualizeResults(config, subject, query)


        print("Compressing results")
        #delete log folder
        shutil.rmtree("ProDA/log")
        #compress the results
        shutil.make_archive("ProDA", "zip", "ProDA")
        print("Finished ProDA")
    #catch and print the error if one is thrown
    except Exception as ex:
        print("\t\t\t\t\033[91mError while running ProDA\033[0m")
        print(ex)
