import argparse
import yaml
import os
import subprocess
from Bio import SeqIO
import csv
import traceback
import multiprocessing as mp
from functools import partial
import math
import pandas as pd
import shutil
import re

manager = mp.Manager()
index = manager.Value("i", -1)
lock = manager.Lock()

def readConfigFile(config_path):
    config_data = None
    with open(config_path, "r") as config_reader:
        config_data = yaml.safe_load(config_reader)

    return config_data

def buildBLASTDatabase(config, subject):
    db_path = "blast_dbs/" + subject + "/"
    log_path = "log/blast_dbs/" + subject + "/"
    if(not os.path.exists(db_path)):
        os.makedirs(db_path)
    if(not os.path.exists(log_path)):
        os.makedirs(log_path)

    files_process = subprocess.Popen("ls -p " + db_path + " | grep -v /$",
                    stdout=subprocess.PIPE, shell=True)
    db_files = list(filter(None, files_process.communicate()[0].decode("utf-8").split("\n")))
    if(not (subject+".nhr" in db_files and subject+".nin" in db_files and subject+".nsq" in db_files)):
        print("Building database for subject '" + subject + "'")
        db_output = subprocess.Popen("(makeblastdb -in " + config["subjects"][subject] + ""
                    " -dbtype nucl -out " + db_path + subject + ") 2> " + log_path + "makeblastdb.log",
                    stdout=subprocess.PIPE, shell=True)
        res,err = db_output.communicate()
        if(res != None):
            print("Building database for subject '" + subject + "'\tfinished")
        elif(err != None):
            print("\033[91mError while building database for subject '" + subject + "',"
                  " see log file '" + log_path + "'\033[0m")
        else:
            print("Found database files for subject '" + subject + "', skipping")

def runBLASTCommand(config, subject, query):
    blast_path = "blast_results/" + subject + "/"
    log_path = "log/blast_results/" + subject + "/"
    if(not os.path.exists(blast_path)):
        os.makedirs(blast_path)
    if(not os.path.exists(log_path)):
        os.makedirs(log_path)

    db_path = "blast_dbs/" + subject + "/" + subject
    if(not os.path.exists(blast_path + query + ".fna")):
        print("Matcher: subject: " + subject + "\tquery: " + query)
        if(not os.path.exists(blast_path + query + ".hit")):
            print("\tBLAST: evalue: " + str(config["evalue"]) + "\tthreads: " + str(config["threads"]))
            blast_output = subprocess.Popen("(tblastn -query " + config["queries"][query] + ""
                           " -db "+ db_path + " -outfmt '6 qseqid sseqid sstart send evalue'"
                           " -evalue " + str(config["evalue"]) + ""
                           " -num_threads 1"
                           " > " + blast_path + query + ".hit" + ")"
                           " 2> " + log_path + query + ".log",
                           stdout=subprocess.PIPE, shell=True)
            res,err = blast_output.communicate()
            if(res != None):
                print("\tBLAST: finished")
            if(err != None):
                print("\t\033[91mError for BLAST: subject: " + subject + "\tquery: " + query + ","
                      " see log file '" + log_path + "'\033[0m")
        else:
            print("\tFound BLAST file, skipping")

        matchBLASTPositions(blast_path + query + ".hit", config["subjects"][subject],
                            config["left_addendum"], config["right_addendum"], log_path)
        print("Matcher: subject: " + subject + "\tquery: " + query + "\tfinished")
    else:
        print("Found matching sequence file for subject '" + subject + "' and"
              " query '" + query + "', skipping")

def matchBLASTPositions(hit_path, subject, lad, rad, log):
    if(os.stat(hit_path).st_size != 0):
        log_path = log + hit_path.split("/")[-1].split(".hit")[0] + ".log"
        try:
            subject_fasta = SeqIO.index(subject, "fasta")
            fasta_results = []
            with open(hit_path, "r") as hit_reader:
                hit_content = csv.reader(hit_reader, delimiter="\t")
                current_query = None
                current_subject = None
                current_sequence = None
                positions = []
                for row in hit_content:
                    if(row[0] != current_query or row[1] != current_subject):
                        if(len(positions)):
                            start = int(positions[0][0])
                            end = int(positions[-1][1])
                            if(start > end):
                                temp_start = start
                                start = end
                                end = temp_start
                            start -= lad
                            if(start < 0):
                                start = 0
                            end += rad
                            if(end >= len(current_sequence)):
                                end = len(current_sequence) - 1

                            header = ">" + current_query + "::" + current_subject + "::" + str(start) + "-" + str(end)
                            sequence = current_sequence[start:end]
                            fasta_results.append(header + "\n" + sequence)
                            positions = []

                        current_query = row[0]
                        current_subject = row[1]
                        current_sequence = str(subject_fasta[current_subject].seq)

                    positions.append([row[2], row[3]])

                if(len(positions)):
                    start = int(positions[0][0])
                    end = int(positions[-1][1])
                    if(start > end):
                        temp_start = start
                        start = end
                        end = temp_start
                    start -= lad
                    if(start < 0):
                        start = 0
                    end += rad
                    if(end >= len(current_sequence)):
                        end = len(current_sequence) - 1

                    header = ">" + current_query + "::" + current_subject + "::" + str(start) + "-" + str(end)
                    sequence = current_sequence[start:end]
                    fasta_results.append(header + "\n" + sequence)
                    del positions[:]
            with open(hit_path.replace(".hit", ".fna"), "w") as fasta_writer:
                fasta_writer.write("\n".join(fasta_results))
            del fasta_results[:]
        except:
            print("\t\033[91mError while matching sequences, see log file '" + log_path + "'\033[0m")
            with open(log_path, "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))
    else:
        with open(hit_path.replace(".hit", ".fna"), "w") as fasta_writer:
            fasta_writer.write("")

def runExonerateCommand(config, subject, query):
    ex_output = "alignment/exonerate/" + subject + "/"
    log_path = "log/alignment/exonerate/" + subject + "/"
    if(not os.path.exists(ex_output)):
        os.makedirs(ex_output)
    if(not os.path.exists(log_path)):
        os.makedirs(log_path)

    print("Exonerate: subject: " + subject + "\tquery: " + query)
    blast_seq = "blast_results/" + subject + "/" + query + ".fna"
    if(not (os.path.exists(ex_output + query + ".gff") and os.path.exists(ex_output + query + ".cds"))):
        if(not os.path.exists(ex_output + query + ".ryo")):
            if(os.path.exists(blast_seq)):
                if(os.stat(blast_seq).st_size != 0):
                    index.value = -1
                    matches = list(SeqIO.parse(blast_seq, "fasta"))
                    pool = mp.Pool(processes=config["threads"])
                    pool_map = partial(runExonerateMultiprocessing, queries=config["queries"][query],
                                       percent=config["exonerate_percentage"], blosum=config["blosum"],
                                       output=ex_output, log=log_path)
                    ryo_results = pool.map_async(pool_map, matches)
                    pool.close()
                    pool.join()
                    ryo_joined_results = [r.strip() for result in ryo_results.get() for r in result]
                    with open(ex_output + query + ".ryo", "w") as ryo_writer:
                        ryo_writer.write("\n".join(ryo_joined_results))

                    del ryo_joined_results[:]
                    del ryo_results.get()[:]
                    print("\tExonerate: finished")
                else:
                    with open(ex_output + query + ".ryo", "w") as ryo_writer:
                        ryo_writer.write("")
                    with open(ex_output + query + ".gff", "w") as gff_writer:
                        gff_writer.write("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")
                    with open(ex_output + query + ".cds", "w") as cds_writer:
                        cds_writer.write("")
            else:
                print("\tRequired files for Exonerate not found, starting BLAST")
                runBLASTCommand(config, subject, query)
                runExonerateCommand(config, subject, query)
        else:
            print("\tExonerate file found, skipping")
            print("\tExonerate processing")
            parseExonerateResults(ex_output + query + ".ryo")
            print("\tExonerate processing finished")
            print("Exonerate: subject: " + subject + "\tquery: " + query + "\tfinished")
    else:
        print("\tProcessed Exonerate files found, skipping")

def runExonerateMultiprocessing(match, queries, percent, blosum, output, log):
    print("\tTarget: " + match.id + "\tpercentage: " + str(percent) + "\tblosum: " + str(blosum))
    local_index = None
    with lock:
        index.value += 1
        local_index = index.value

    query_fasta = SeqIO.index(queries, "fasta")
    query = match.id.split("::")[0]
    temp_target = output + query + "_" + str(local_index) + "_target.fna"
    with open(temp_target, "w") as target_writer:
        target_writer.write(">" + match.id + "\n" + str(match.seq))

    temp_query = output + query + "_" + str(local_index) + "_query.fna"
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query_fasta[query].id + "\n" + str(query_fasta[query].seq))

    temp_output_ryo = output + query + "_" + str(local_index) + ".ryo"
    ex_output = subprocess.Popen("(exonerate --model protein2genome --targettype dna --querytype protein"
                " --showtargetgff --showalignment no --showvulgar no --ryo '>%ti\n%tcs'"
                " --refine region --proteinsubmat blosum/" + str(config["blosum"]) + ".txt"
                " --percent " + str(config["exonerate_percentage"]) + ""
                " --query " + temp_query + " --target " + temp_target + ""
                " > " + temp_output_ryo + ")"
                " 2> " + log + query + ".log",
                stdout=subprocess.PIPE, shell=True)
    res,err = ex_output.communicate()
    if(err != None):
        print("\t\033[91mError while running Exonerate, see log file '" + log + query + ".log" + "'\033[0m")

    temp_output_ryo2 = output + query + "_" + str(local_index) + ".ryo2"
    os.system("(tail -n +4 " + temp_output_ryo + " | head -n -1) > " + temp_output_ryo2)
    ryo_results = None
    with open(temp_output_ryo2, "r") as output_reader:
        ryo_results = output_reader.readlines()

    os.remove(temp_target)
    os.remove(temp_query)
    os.remove(temp_output_ryo)
    os.remove(temp_output_ryo2)
    print("\tTarget: " + match.id + "\tfinished")
    return ryo_results

def parseExonerateResults(ryo_results):
    if(os.stat(ryo_results).st_size  != 0):
        gff = []
        cds = []
        seq_found = False
        positions = []
        with open(ryo_results, "r") as ryo_reader:
            lines = ryo_reader.readlines()
            hit_id = 0
            for line in lines:
                if(line):
                    if(line.startswith("#")):
                        seq_found = False
                    elif(line.startswith(">")):
                        header = line.strip().split("::")[0] + "::" + line_splitted[0].split("::")[1] + "::" + str(hit_id)
                        cds.append(header)
                        hit_id += 1
                        seq_found = True
                    elif(not seq_found):
                        line_splitted = line.strip().split("\t")
                        query = line_splitted[0].split("::")[0]
                        target = line_splitted[0].split("::")[1]
                        start = int(line_splitted[3]) + int(line_splitted[0].split("::")[2].split("-")[0]) - 1
                        end = int(line_splitted[4]) + int(line_splitted[0].split("::")[2].split("-")[0])
                        if(line_splitted[2] == "gene"):
                            if(len(positions)):
                                if(int(positions[0].split("\t")[5]) > int(positions[-1].split("\t")[5])):
                                    gff.append("\n".join(list(reversed(positions))))
                                else:
                                    gff.append("\n".join(positions))
                            positions = []
                            gff.append("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")
                            gff.append(target + "\t" + query + "\texonerate:aln\t" + line_splitted[2] + "\t" + str(start) + "\t" + str(end) + "\t" + line_splitted[6] + "\t" + line_splitted[5] + "\t" + str(hit_id))
                        elif(line_splitted[2] == "cds"):
                            positions.append(target + "\t" + query + "\texonerate:aln\t" + line_splitted[2] + "\t" + str(start) + "\t" + str(end) + "\t" + line_splitted[6] + "\t" + line_splitted[5] + "\t" + str(hit_id))
                    elif(seq_found):
                        cds.append(re.sub(r'[^A-Za-z*]', '', line.strip()))
            if(len(positions)):
                if(int(positions[0].split("\t")[5]) > int(positions[-1].split("\t")[5])):
                    gff.append("\n".join(list(reversed(positions))))
                else:
                    gff.append("\n".join(positions))

        with open(ryo_results.replace(".ryo", ".gff"), "w") as gff_writer:
            gff_writer.write("\n".join(gff))
        with open(ryo_results.replace(".ryo", ".cds"), "w") as cds_writer:
            cds_writer.write("\n".join(cds))

        del gff[:]
        del cds[:]
        del positions[:]
    else:
        with open(ryo_results.replace(".ryo", ".gff"), "w") as gff_writer:
            gff_writer.write("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")
        with open(ryo_results.replace(".ryo", ".cds"), "w") as cds_writer:
            cds_writer.write("")

def runSpalnCommand(config, subject, query):
    sp_output = "alignment/spaln/" + subject + "/"
    log_path = "log/alignment/spaln/" + subject + "/"
    if(not os.path.exists(sp_output)):
        os.makedirs(sp_output)
    if(not os.path.exists(log_path)):
        os.makedirs(log_path)

    print("Spaln: subject: " + subject + "\tquery: " + query)
    blast_seq = "blast_results/" + subject + "/" + query + ".fna"
    if(not (os.path.exists(sp_output + query + ".gff") and os.path.exists(sp_output + query + ".cds"))):
        if(not os.path.exists(sp_output + query + ".sp")):
            if(os.path.exists(blast_seq)):
                if(os.stat(blast_seq).st_size != 0):
                    index.value = -1
                    matches = list(SeqIO.parse(blast_seq, "fasta"))
                    pool = mp.Pool(processes=config["threads"])
                    pool_map = partial(runSpalnMultiprocessing, queries=config["queries"][query],
                                       pam=config["pam"], output=sp_output, log=log_path)
                    sp_results = pool.map_async(pool_map, matches)
                    pool.close()
                    pool.join()
                    sp_joined_results = [r.strip() for result in sp_results.get() for r in result]
                    with open(sp_output + query + ".sp", "w") as sp_writer:
                        sp_writer.write("\n".join(sp_joined_results))

                    del sp_joined_results[:]
                    del sp_results.get()[:]
                    print("\tSpaln: finished")
                    print("\tSpaln processing")
                    parseSpalnResults(sp_output + query + ".sp")
                    print("\tSpaln processing finished")
                    print("Spaln: subject: " + subject + "\tquery: " + query + "\tfinished")
                else:
                    with open(sp_output + query + ".sp", "w") as sp_writer:
                        sp_writer.write("")
                    with open(sp_output + query + ".gff", "w") as sp_writer:
                        sp_writer.write("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")
                    with open(sp_output + query + ".cds", "w") as sp_writer:
                        sp_writer.write("")
            else:
                print("\tRequired files for Spaln not found, starting BLAST")
                runBLASTCommand(config, subject, query)
                runSpalnCommand(config, subject, query)
        else:
            print("\tSpaln file found, skipping")
            print("\tSpaln processing")
            parseSpalnResults(sp_output + query + ".sp")
            print("\tSpaln processing finished")
            print("Spaln: subject: " + subject + "\tquery: " + query + "\tfinished")
    else:
        print("\tProcessed Spaln files found, skipping")

def runSpalnMultiprocessing(match, queries, pam, output, log):
    print("\tTarget: " + match.id + "\tpam: " + str(pam))
    local_index = None
    with lock:
        index.value += 1
        local_index = index.value

    query_fasta = SeqIO.index(queries, "fasta")
    query = match.id.split("::")[0]
    temp_target = output + query + "_" + str(local_index) + "_target.fna"
    with open(temp_target, "w") as target_writer:
        target_writer.write(">" + match.id + "\n" + str(match.seq))

    temp_query = output + query + "_" + str(local_index) + "_query.fna"
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query_fasta[query].id + "\n" + str(query_fasta[query].seq))

    temp_output = output + query + "_" + str(local_index) + ".sp"
    spaln_output = subprocess.Popen("($HOME/spaln2.4.0/bin/spaln -M4 -N1 -Q3 -S3 -yp" + str(pam) + " -yq" + str(pam) + ""
                   " -O6 -o" + temp_output + " " + temp_target + " " + temp_query + ")"
                   " 2> " + log + query + ".log",
                   stdout=subprocess.PIPE, shell=True)
    res,err = spaln_output.communicate()
    if(err != None):
        print("\t\033[91mError while running Spaln, see log file '" + log + query + ".log'\033[0m")
    sp_results = None
    with open(temp_output, "r") as output_reader:
        sp_results = output_reader.readlines()

    os.remove(temp_target)
    os.remove(temp_query)
    os.remove(temp_output)
    print("\tTarget: " + match.id + "\tfinished")
    return sp_results

def parseSpalnResults(sp_results):
    if(os.stat(sp_results).st_size != 0):
        gff = []
        cds = []
        with open(sp_results, "r") as sp_reader:
            lines = sp_reader.readlines()
            target = None
            query = None
            orientation = None
            block_start = None
            pos = []
            hit_id = 0
            for line in lines:
                if(line and not line.startswith(";M")):
                    if(line.startswith(">")):
                        if(len(pos)):
                            for position in pos:
                                    for start_end_pos in position:
                                        start = int(start_end_pos.split("..")[0].strip()) + block_start - 1
                                        end = int(start_end_pos.split("..")[1].strip()) + block_start
                                        gff.append(target + "\t" + query +"\tspaln:aln\tcds\t" + str(start) + "\t" + str(end) + "\t" + orientation + "\t.\t" + str(hit_id))
                            pos = []
                            hit_id += 1

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
                        if(start > end):
                            temp_start = start
                            start = end
                            end = temp_start
                        start -= 1

                        gff.append(target + "\t" + query + "\tspaln:aln\tgene\t" + str(start) + "\t" + str(end) + "\t" + orientation + "\t" + score + "\t" + str(hit_id))
                        cds.append(header)
                    elif(line.startswith(";C")):
                        if("join" in line):
                            pos.append(list(filter(None, line.strip().split("join(")[1].replace(")", "").split(","))))
                        else:
                            pos.append(list(filter(None, line.strip().split(";C ")[1].replace(")", "").split(","))))
                    else:
                        cds.append(re.sub(r'[^A-Za-z*]', '', line.strip()))

            if(len(pos)):
                for position in pos:
                        for start_end_pos in position:
                            start = int(start_end_pos.split("..")[0].strip()) + block_start - 1
                            end = int(start_end_pos.split("..")[1].strip()) + block_start
                            gff.append(target + "\t" + query + "\tspaln:aln\tcds\t" + str(start) + "\t" + str(end) + "\t" + orientation + "\t.\t" + str(hit_id))
                del pos[:]

        with open(sp_results.replace(".sp", ".gff"), "w") as gff_writer:
            gff_writer.write("\n".join(gff))
        with open(sp_results.replace(".sp", ".cds"), "w") as cds_writer:
            cds_writer.write("\n".join(cds))

        del gff[:]
        del cds[:]
    else:
        with open(sp_results.replace(".sp", ".gff"), "w") as gff_writer:
            gff_writer.write("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")
        with open(sp_results.replace(".sp", ".cds"), "w") as cds_writer:
            cds_writer.write("")

def evaluateAlignment(config, algorithm, subject, query):
    output = "evaluation/" + algorithm + "/" + subject + "/"
    log_path = "log/evaluation/" + algorithm + "/" + subject + "/"
    if(not os.path.exists(output)):
        os.makedirs(output)
    if(not os.path.exists(log_path)):
        os.makedirs(log_path)

    print("Evaluation: subject: " + subject + "\tquery: " + query)
    hits_path = "alignment/" + algorithm + "/" + subject + "/" + query + ".cds"
    gff = "alignment/" + algorithm + "/" + subject + "/" + query + ".gff"
    if(not os.path.exists(output + query + ".stretcher")):
        if(os.path.exists(hits_path) and os.path.exists(gff)):
            if(os.stat(hits_path).st_size != 0):
                index.value = -1
                hits = list(SeqIO.parse(hits_path, "fasta"))
                pool = mp.Pool(processes=config["threads"])
                pool_map = partial(runStretcherMultiprocessing, queries=config["queries"][query],
                                   distance=config["hssp_distance"], stop=config["stop_at_stop"],
                                   len_cutoff=config["len_cutoff"], id_thres=config["identity_threshold"],
                                   m_fil=config["filter_met"], output=output, log=log_path)
                stretcher_results = pool.map_async(pool_map, hits)
                pool.close()
                pool.join()
                with open(output + query + ".stretcher", "w") as eval_writer:
                    eval_writer.write("#id\ttarget\tquery\tidentity\tsimilarity\taln_length\thssp_identity\n"
                                      "" + "\n".join(filter(None, stretcher_results.get())))
                print("Evaluation: subject: " + subject + "\tquery: " + query + "\tfinished")
            else:
                with open(output + query + ".stretcher", "w") as eval_writer:
                    eval_writer.write("#id\ttarget\tquery\tidentity\tsimilarity\taln_length\thssp_identity\n")
        else:
            if(algorithm == "exonerate"):
                print("\tRequired files for evaluation not found, starting Exonerate")
                runExonerateCommand(config, subject, query)
            else:
                print("\tRequired files for evaluation not found, starting Spaln")
                runSpalnCommand(config, subject, query)
            evaluateAlignment(config, algorithm, subject, query)
    else:
        print("\tEvaluation file found, skipping")

def runStretcherMultiprocessing(hit, queries, distance, stop, len_cutoff, id_thres, m_fil, output, log):
    print("\tTarget: " + hit.id + "\tstop_at_stop: " + str(stop) + ""
          "\tlength_cutoff_percentage: " + str(len_cutoff) + "\tidentity_threshold: " + str(id_thres) + ""
          "\tHSSP_distance: " + str(distance) + "\tmet_filter: " + str(m_fil))
    local_index = None
    with lock:
        index.value += 1
        local_index = index.value

    query_fasta = SeqIO.index(queries, "fasta")
    query = hit.id.split("::")[0]
    try:
        subject = hit.id.split("::")[1]
        id = hit.id.split("::")[2].strip()
        temp_target = output + query + "_" + str(local_index) + "_target.faa"

        translated_sequence = None
        if(stop > 0):
            translated_sequence = str(trimSequence(hit.seq).translate(to_stop=True))
        else:
            translated_sequence = str(trimSequence(hit.seq).translate())

        with open(temp_target, "w") as target_writer:
            target_writer.write(">" + hit.id + "\n" + translated_sequence)

        temp_query = output + query + "_" + str(local_index) + "_query.faa"
        query_seq = str(query_fasta[query].seq)
        with open(temp_query, "w") as query_writer:
            query_writer.write(">" + query_fasta[query].id + "\n" + query_seq)

        st_result = []
        if(m_fil <= 0 or translated_sequence.startswith("M")):
            temp_output = output + query + "_" + str(local_index) + ".stretcher"
            os.system("(stretcher -asequence " + temp_query + " -sprotein1"
                      " -bsequence " + temp_target + " -sprotein2 -auto -stdout > " + temp_output + ")"
                      " 2> " + log + query + ".log")
            identity = None
            al_length = None
            identity = None
            similarity = None
            with open(temp_output, "r") as output_reader:
                    content = output_reader.readlines()
                    for line in content:
                        if(line):
                            if(line.startswith("# Length")):
                                al_length = int(line.split(" ")[2].strip())
                            if(line.startswith("# Identity")):
                                identity = float(line.split("(")[1][:-3].strip())
                            if(line.startswith("# Similarity")):
                                similarity = float(line.split("(")[1][:-3].strip())
                                break
            hssp_identity = calculateHSSPIdentity(al_length, distance)
            if(similarity >= identity and identity >= hssp_identity and identity >= id_thres):
                st_result = [id, subject, query, str(identity), str(similarity), str(al_length), str(hssp_identity)]
            os.remove(temp_output)

        os.remove(temp_target)
        os.remove(temp_query)
        print("\tTarget: " + hit.id + "\tfinished")
    except:
        print("\t\033[91mError while evaluating, see log file '" + log + query + ".log" + "'\033[0m")
        with open(log + query + ".log", "w") as log_writer:
            log_writer.write(str(traceback.format_exc()))
    return("\t".join(st_result))

def trimSequence(seq):
    return seq[:len(seq)-len(seq)%3]

def calculateHSSPIdentity(length, distance):
    if(length < 12):
        return 100
    elif(length < 418):
        return distance + 480 * length**(-0.32 * (1 + math.exp(-length / 1000)))
    else:
        return distance + 19.5

def mergeResults(config, subject, query):
        subject_fasta = SeqIO.index(config["subjects"][subject], "fasta")
        print("Merging: subject: " + subject + "\tquery: " + query + "\toverlap percent: " + str(config["overlap_percentage"]))
        try:
            output = "results/" + subject + "/" + query + "/"
            log_path = "log/results/" + subject + "/" + query + "/"
            if(not os.path.exists(output)):
                os.makedirs(output)
            if(not os.path.exists(log_path)):
                os.makedirs(log_path)
            gff_list = ["#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid"]
            cds_list = []
            ex_gff = "alignment/exonerate/" + subject + "/" + query + ".gff"
            ex_stretcher = "evaluation/exonerate/" + subject + "/" + query + ".stretcher"
            ex_id_list = []
            sp_gff = "alignment/spaln/" + subject + "/" + query + ".gff"
            sp_stretcher = "evaluation/spaln/" + subject + "/" + query + ".stretcher"
            sp_id_list = []
            hit_id = 0
            if(not (os.path.exists(output + query + ".gff") or os.path.exists(output + query + ".nuc") or os.path.exists(output + query + ".cds") or os.path.exists(output + query + ".pep"))):
                if(os.path.exists(ex_gff) and os.path.exists(ex_stretcher) and os.path.exists(sp_gff) and os.path.exists(sp_stretcher)):
                    sp_sequences = SeqIO.index("alignment/spaln/" + subject + "/" + query + ".cds", "fasta")
                    sp_gff_df = pd.read_csv(sp_gff, sep="\t", header=0).astype(str)
                    sp_eval_df = pd.read_csv(sp_stretcher, sep="\t", header=0).astype(str)
                    ex_sequences = SeqIO.index("alignment/exonerate/" + subject + "/" + query + ".cds", "fasta")
                    ex_gff_df = pd.read_csv(ex_gff, sep="\t", header=0).astype(str)
                    ex_eval_df = pd.read_csv(ex_stretcher, sep="\t", header=0).astype(str)
                    for sp_index,sp_row in sp_eval_df.iterrows():
                        sp_header = sp_row["query"] + "::" + sp_row["target"] + "::" + sp_row["#id"]
                        sp_seq = str(sp_sequences[sp_header].seq)
                        ex_gene_rows = ex_gff_df.loc[(ex_gff_df["type"] == "gene") & (ex_gff_df["#subject"] == sp_row["target"]) & (ex_gff_df["query"] == sp_row["query"])]
                        ex_eval_rows = ex_eval_df.loc[(ex_eval_df["target"] == sp_row["target"]) & (ex_eval_df["query"] == sp_row["query"])]
                        sp_gene_rows = sp_gff_df.loc[(sp_gff_df["type"] == "gene") & (sp_gff_df["#subject"] == sp_row["target"]) & (sp_gff_df["query"] == sp_row["query"])]
                        sp_id_row = sp_gene_rows.loc[sp_gene_rows["id"] == sp_row["#id"]].iloc[0]
                        for ex_index,ex_row in ex_eval_rows.iterrows():
                            if(not ex_row["#id"] in ex_id_list):
                                ex_header = ex_row["query"] + "::" + ex_row["target"] + "::" + ex_row["#id"]
                                ex_seq = str(ex_sequences[ex_header].seq)
                                ex_id_row = ex_gene_rows.loc[ex_gene_rows["id"] == ex_row["#id"]].iloc[0]
                                if(calculateOverlapPercentage(int(sp_id_row["start"]), int(sp_id_row["end"]), int(ex_id_row["start"]), int(ex_id_row["end"])) >= config["overlap_percentage"]):
                                    sp_id_list.append(sp_row["#id"])
                                    ex_id_list.append(ex_row["#id"])
                                    if((sp_seq.startswith("ATG") and ex_seq.startswith("ATG"))
                                       or not (sp_seq.startswith("ATG") or ex_seq.startswith("ATG"))):
                                       if(float(sp_id_row["score"]) >= float(ex_id_row["score"])):
                                           id_rows = sp_gff_df.loc[sp_gff_df["id"] == str(sp_row["#id"])]
                                           id_rows["id"] = str(hit_id)
                                           header = "##" +  sp_row["target"] + "\t" + sp_row["query"] + "\t" + str(hit_id)
                                           gff_list.append(header + "\n" + id_rows.to_string(header=False, index=False))
                                           cds_list.append(">" + sp_row["query"] + "::" + sp_row["target"] + "::" + str(hit_id) + "\n" + sp_seq)
                                           hit_id += 1
                                       else:
                                           id_rows = ex_gff_df.loc[ex_gff_df["id"] == ex_row["#id"]]
                                           id_rows["id"] = str(hit_id)
                                           header = "##" +  ex_row["target"] + "\t" + ex_row["query"] + "\t" + str(hit_id)
                                           gff_list.append(header + "\n" + id_rows.to_string(header=False, index=False).strip())
                                           cds_list.append(">" + ex_row["query"] + "::" + ex_row["target"] + "::" + str(hit_id) + "\n" + ex_seq)
                                           hit_id += 1
                                    elif(sp_seq.startswith("ATG")):
                                        id_rows = sp_gff_df.loc[sp_gff_df["id"] == sp_row["#id"]]
                                        id_rows["id"] = hit_id
                                        header = "##" +  sp_row["target"] + "\t" + sp_row["query"] + "\t" + str(hit_id)
                                        gff_list.append(header + "\n" + id_rows.to_string(header=False, index=False))
                                        cds_list.append(">" + sp_row["query"] + "::" + sp_row["target"] + "::" + str(hit_id) + "\n" + sp_seq)
                                        hit_id += 1
                                    else:
                                        id_rows = ex_gff_df.loc[ex_gff_df["id"] == ex_row["#id"]]
                                        id_rows["id"] = str(hit_id)
                                        header = "##" +  ex_row["target"] + "\t" + ex_row["query"] + "\t" + str(hit_id)
                                        gff_list.append(header + "\n" + id_rows.to_string(header=False, index=False))
                                        cds_list.append(">" + ex_row["query"] + "::" + ex_row["target"] + "::" + str(hit_id) + "\n" + ex_seq)
                                        hit_id += 1
                                    break

                    for sp_index,sp_row in sp_eval_df.iterrows():
                        if(not sp_row["#id"] in sp_id_list):
                            id_rows = sp_gff_df.loc[sp_gff_df["id"] == sp_row["#id"]]
                            id_rows["id"] = str(hit_id)
                            header = "##" +  sp_row["target"] + "\t" + sp_row["query"] + "\t" + str(hit_id)
                            sp_header = sp_row["query"] + "::" + sp_row["target"] + "::" + sp_row["#id"]
                            gff_list.append(header + "\n" + id_rows.to_string(header=False, index=False))
                            cds_list.append(">" + sp_row["query"] + "::" + sp_row["target"] + "::" + str(hit_id) + "\n" + str(sp_sequences[sp_header].seq))
                            hit_id += 1

                    for ex_index,ex_row in ex_eval_df.iterrows():
                        if(not ex_row["#id"] in ex_id_list):
                            id_rows = ex_gff_df.loc[ex_gff_df["id"] == ex_row["#id"]]
                            id_rows["id"] = str(hit_id)
                            header = "##" +  ex_row["target"] + "\t" + ex_row["query"] + "\t" + str(hit_id)
                            ex_header = ex_row["query"] + "::" + ex_row["target"] + "::" + ex_row["#id"]
                            gff_list.append(header + "\n" + id_rows.to_string(header=False, index=False))
                            cds_list.append(">" + ex_row["query"] + "::" + ex_row["target"] + "::" + str(hit_id) + "\n" + str(ex_sequences[ex_header].seq))
                            hit_id += 1

                    with open(output + query + ".gff", "w") as gff_writer:
                        gff_writer.write("\n".join(gff_list))

                    nucl = []
                    with open(output + query + ".gff", "r") as gff_reader:
                        content = csv.reader(gff_reader, delimiter="\t")
                        target = None
                        protein = None
                        id = None
                        for row in content:
                            if(row[0].startswith("##")):
                                target = row[0][2:]
                                protein = row[1]
                                id = row[2]
                            elif(len(row) == 1 and row[0].split(" ")[7] == "gene"):
                                row_filter = list(filter(None, row[0].split(" ")))
                                start = int(row_filter[4])
                                end = int(row_filter[5])
                                nucl.append(">" + protein + "::" + target + "::" + id + "\n" + str(subject_fasta[target].seq)[start:end])
                    with open(output + query + ".nuc", "w") as cds_writer:
                        cds_writer.write("\n".join(nucl))

                    with open(output + query + ".cds", "w") as cds_writer:
                        cds_writer.write("\n".join(cds_list))

                    translated_cds = []
                    for fasta in SeqIO.parse(output + query + ".cds", "fasta"):
                        translated_cds.append(">" + fasta.id + "\n" + str(fasta.seq.translate()))
                    with open(output + query + ".pep", "w") as cds_writer:
                        cds_writer.write("\n".join(translated_cds))

                    del gff_list[:]
                    del nucl[:]
                    del cds_list[:]
                    del translated_cds[:]
                    del ex_id_list[:]
                    del sp_id_list[:]
                    print("Merging: subject: " + subject + "\tquery: " + query + "\tfinished")
                else:
                    if(not os.path.exists(ex_gff)):
                        print("\tRequired Exonerate files for Merger not found, starting Exonerate")
                        runExonerateCommand(config, subject, query)
                    if(not os.path.exists(sp_gff)):
                        print("\tRequired Spaln files for Merger not found, starting Spaln")
                        runSpalnCommand(config, subject, query)
                    if(not os.path.exists(ex_stretcher)):
                        print("\tRequired Ex_Eval files for Merger not found, starting Exonerate evaluation")
                        evaluateAlignment(config, "exonerate", subject, query)
                    if(not os.path.exists(sp_stretcher)):
                        print("\tRequired Sp_Eval files for Merger not found, starting Spaln evaluation")
                        evaluateAlignment(config, "spaln", subject, query)
                    mergeResults(config, subject, query)
            else:
                print("\tResult files found, skipping")
        except:
            print("\t\033[91mError while merging results, see log file '" + log_path + query + ".log'\033[0m")
            with open(log_path + query + ".log", "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))

def calculateOverlapPercentage(start1, end1, start2, end2):
    if(end1 < start2 or end2 < start1):
        return 0
    else:
        return ((min(end1, end2) - max(start1, start2)) / min(end1-start1, end2-start2)) * 100

def visualizeResults(config, subject, query):
    print("Visualizing results, subject: " + subject + "\tquery: " + query + "\tprop file: " + config["properties"])
    protein = "results/" + subject + "/" + query + "/" + query + ".pep"
    output = "results/" + subject + "/" + query + "/alignment/"
    if(not os.path.exists(output)):
        if(os.path.exists(protein)):
            index.value = -1
            os.makedirs(output)
            log_path = "log/results/" + subject + "/" + query + "/alignment/"
            if(not os.path.exists(log_path)):
                os.makedirs(log_path)
            protein_fasta = list(SeqIO.parse(protein, "fasta"))
            pool = mp.Pool(processes=config["threads"])
            pool_map = partial(runJalViewMultiprocessing, queries=config["queries"][query],
                               properties=config["properties"], output=output, log=log_path)
            pool.map_async(pool_map, protein_fasta)
            pool.close()
            pool.join()
            print("Visualizing: subject: " + subject + "\tquery: " + query + "\tfinished")
        else:
            print("\tRequired files for Visualizer not found, starting Merger")
            mergeResults(config, subject, query)
            visualizeResults(config, subject, query)
    else:
        print("\tVisualizing folder found, skipping")

def runJalViewMultiprocessing(hit, queries, properties, output, log):
    print("\tTarget: " + hit.id)
    local_index = None
    with lock:
        index.value += 1
        local_index = index.value

    query_fasta = SeqIO.index(queries, "fasta")
    query = hit.id.split("::")[0]
    id = hit.id.split("::")[2]
    temp_target = output + query + "_" + str(local_index) + "_target.fna"
    with open(temp_target, "w") as target_writer:
        target_writer.write(">" + hit.id + "\n" + str(hit.seq))

    temp_query = output + query + "_" + str(local_index) + "_query.fna"
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query_fasta[query].id + "\n" + str(query_fasta[query].seq))

    clustal = output + query + "_" + id + ".clustal"
    os.system("(stretcher -asequence " + temp_query + " -sprotein1 -bsequence " + temp_target + ""
              " -sprotein1 -auto -aformat clustal -stdout > " + clustal + ") 2> " + log + query + "_" + id + ".cl_log")
    png = output + query + "_" + id + ".png"
    jalview_output = subprocess.Popen("(jalview -nodisplay -props " + properties + " -colour clustal"
                     " -open " + clustal + " -png " + png + ") 2> " + log + query + "_" + id + ".jv_log",
                     stdout=subprocess.PIPE, shell=True)
    res,err = jalview_output.communicate()
    if(err != None):
        print("\t\033[91mError while visualizing results, see log file '" + log + query + "_" + id + ".jv_log'\033[0m")

    os.remove(temp_target)
    os.remove(temp_query)
    print("\tTarget: " + hit.id + "\tfinished")


parser = argparse.ArgumentParser(description="ProDA (Protein-DNA-Aligner): Search for protein sequences in DNA sequences with the Exonerate and Spaln alignment tools and analyze the results.")
parser.add_argument("-c", "--config", help="Path to config.yaml file", default="config.yaml")
args = parser.parse_args()
if(not os.path.exists(args.config)):
    print("\033[91mError: Unable to find/access config file named '" + args.config + "'\033[0m.")
else:
    config = readConfigFile(args.config)
    for subject in config["subjects"]:
        buildBLASTDatabase(config, subject)
        for query in config["queries"]:
            visualizeResults(config, subject, query)
    print("Finished ProDA")
