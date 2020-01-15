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

def readConfigFile(config_path):
    config_data = None
    with open(config_path, "r") as config_reader:
        config_data = yaml.safe_load(config_reader)

    return config_data

def buildBLASTDatabase(config):
    for subject in config["subjects"]:
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
                print("Error while building database for subject '" + subject + "',"
                      " see log file '" + log_path + "'")
        else:
            print("Found database files for subject '" + subject + "', skipping")

def runBLASTCommand(config):
    for subject in config["subjects"]:
        blast_path = "blast_results/" + subject + "/"
        log_path = "log/blast_results/" + subject + "/"
        if(not os.path.exists(blast_path)):
            os.makedirs(blast_path)
        if(not os.path.exists(log_path)):
            os.makedirs(log_path)

        db_path = "blast_dbs/" + subject + "/" + subject
        for query in config["queries"]:
            if(not os.path.exists(blast_path + query + ".fna")):
                print("Matcher: subject: " + subject + "\tquery: " + query)
                if(not os.path.exists(blast_path + query + ".hit")):
                    print("\tBLAST: evalue: " + str(config["evalue"] + "\tthreads: " + str(config["threads"])))
                    blast_output = subprocess.Popen("(tblastn -query " + config["queries"][query] + ""
                                   " -db "+ db_path + " -outfmt '6 qseqid sseqid sstart send evalue'"
                                   " -evalue " + str(config["evalue"]) + ""
                                   " -num_threads " + str(config["threads"]) + ""
                                   " > " + blast_path + query + ".hit" + ")"
                                   " 2> " + log_path + query + ".log",
                                   stdout=subprocess.PIPE, shell=True)
                    res,err = blast_output.communicate()
                    if(res != None):
                        print("\tBLAST: finished")
                    if(err != None):
                        print("\tError for BLAST: subject: " + subject + "\tquery: " + query + ","
                              " see log file '" + log_path + "'")
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
            print("Error while matching sequences, see log file '" + log_path + "'")
            with open(log_path, "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))

def runExonerateCommand(config):
    for subject in config["subjects"]:
        ex_output = "alignment/exonerate/" + subject + "/"
        log_path = "log/alignment/exonerate/" + subject + "/"
        if(not os.path.exists(ex_output)):
            os.makedirs(ex_output)
        if(not os.path.exists(log_path)):
            os.makedirs(log_path)

        for query in config["queries"]:
            print("Exonerate: subject: " + subject + "\tquery: " + query)
            blast_seq = "blast_results/" + subject + "/" + query + ".fna"
            if(os.path.exists(blast_seq)):
                if(os.stat(blast_seq).st_size != 0):
                    if(not (os.path.exists(ex_output + query + ".gff") and os.path.exists(ex_output + query + ".cds"))):
                        if(not os.path.exists(ex_output + query + ".ryo")):
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
                            print("Exonerate: finished")
                        else:
                            print("Exonerate file found, skipping")

                        print("Exonerate processing")
                        parseExonerateResults(ex_output + query + ".ryo")
                        print("Exonerate processing finished")
                    else:
                        print("Processed Exeronate files found, skipping")

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
        print("Error while running Exonerate, see log file '" + log + query + ".log" + "'")

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
    if(os.stat(ryo_results) != 0):
        gff = []
        cds = []
        seq_found = False
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
                            gff.append("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")
                            gff.append(target + "\t" + query + "\texonerate:aln\t" + line_splitted[2] + "\t" + str(start) + "\t" + str(end) + "\t" + line_splitted[6] + "\t" + line_splitted[5] + "\t" + str(hit_id))
                        elif(line_splitted[2] == "cds"):
                            gff.append(target + "\t" + query + "\texonerate:aln\t" + line_splitted[2] + "\t" + str(start) + "\t" + str(end) + "\t" + line_splitted[6] + "\t" + line_splitted[5] + "\t" + str(hit_id))
                    elif(seq_found):
                        cds.append(line.strip())

        with open(ryo_results.replace(".ryo", ".gff"), "w") as gff_writer:
            gff_writer.write("\n".join(gff))

        with open(ryo_results.replace(".ryo", ".cds"), "w") as cds_writer:
            cds_writer.write("\n".join(cds))

        del gff[:]
        del cds[:]

def runSpalnCommand(config):
    for subject in config["subjects"]:
        sp_output = "alignment/spaln/" + subject + "/"
        log_path = "log/alignment/spaln/" + subject + "/"
        if(not os.path.exists(sp_output)):
            os.makedirs(sp_output)
        if(not os.path.exists(log_path)):
            os.makedirs(log_path)

        for query in config["queries"]:
            print("Spaln: subject: " + subject + "\tquery: " + query)
            blast_seq = "blast_results/" + subject + "/" + query + ".fna"
            if(os.path.exists(blast_seq)):
                if(os.stat(blast_seq).st_size != 0):
                    if(not (os.path.exists(sp_output + query + ".gff") and os.path.exists(sp_output + query + ".cds"))):
                        if(not os.path.exists(sp_output + query + ".sp")):
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
                            print("Spaln: finished")
                        else:
                            print("Spaln file found, skipping")

                        print("Spaln processing")
                        parseSpalnResults(sp_output + query + ".sp")
                        print("Spaln processing finished")
                    else:
                        print("Processed Spaln files found, skipping")

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
    spaln_output = subprocess.Popen("(spaln -M -N1 -Q3 -S3 -yp" + str(pam) + " -yq" + str(pam) + ""
                   " -O6 -o" + temp_output + " " + temp_target + " " + temp_query + ")"
                   " 2> " + log + query + ".log",
                   stdout=subprocess.PIPE, shell=True)
    res,err = spaln_output.communicate()
    if(err != None):
        print("Error while running Spaln, see log file '" + log + query + ".log" + "'")

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
                        start = int(line_splitted[6]) + block_start - 1
                        end = int(line_splitted[8]) + block_start
                        orientation = line_splitted[2]
                        score = line_splitted[-1]
                        gff.append("#subject\tquery\talgorithm\ttype\tstart\tend\torientation\tscore\tid")
                        if(start > end):
                            temp_start = start
                            start = end
                            end = temp_start

                        gff.append(target + "\t" + query + "\tspaln:aln\tgene\t" + str(start) + "\t" + str(end) + "\t" + orientation + "\t" + score + "\t" + str(hit_id))
                        cds.append(header)
                    elif(line.startswith(";C")):
                        if("join" in line):
                            pos.append(list(filter(None, line.strip().split("join(")[1].replace(")", "").split(","))))
                        else:
                            pos.append(list(filter(None, line.strip().split(";C ")[1].replace(")", "").split(","))))
                    else:
                        cds.append(line.strip())

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

def evaluateAlignment(config, algorithm):
    for subject in config["subjects"]:
        output = "evaluation/" + algorithm + "/" + subject + "/"
        log_path = "log/evaluation/" + algorithm + "/" + subject + "/"
        if(not os.path.exists(output)):
            os.makedirs(output)
        if(not os.path.exists(log_path)):
            os.makedirs(log_path)

        for query in config["queries"]:
            print("Evaluation: subject: " + subject + "\tquery: " + query)
            hits_path = "alignment/" + algorithm + "/" + subject + "/" + query + ".cds"
            gff = "alignment/" + algorithm + "/" + subject + "/" + query + ".gff"
            if(os.path.exists(hits_path) and os.path.exists(gff)):
                if(not os.path.exists(output + query + ".stretcher")):
                    if(os.stat(hits_path).st_size != 0):
                        hits = list(SeqIO.parse(hits_path, "fasta"))
                        pool = mp.Pool(processes=config["threads"])
                        pool_map = partial(runStretcherMultiprocessing, queries=config["queries"][query],
                                           distance=config["hssp_distance"], use_sim=config["similarity"], output=output, log=log_path)
                        stretcher_results = pool.map_async(pool_map, hits)
                        pool.close()
                        pool.join()
                        with open(output + query + ".stretcher", "w") as eval_writer:
                            eval_writer.write("#id\ttarget\tquery\tidentity\tsimilarity\taln_length\thssp_identity\thssp_similarity\n" + "\n".join(stretcher_results.get()))
                        print("Evaluation: finished")
                else:
                    print("Evaluation file found, skipping")

def runStretcherMultiprocessing(hit, queries, distance, use_sim, output, log):
    print("\tTarget: " + hit.id + "\tuse id: " + str(use_sim))
    local_index = None
    with lock:
        index.value += 1
        local_index = index.value

    query_fasta = SeqIO.index(queries, "fasta")
    query = hit.id.split("::")[0]
    subject = hit.id.split("::")[1]
    id = hit.id.split("::")[2].strip()
    temp_target = output + query + "_" + str(local_index) + "_target.faa"
    with open(temp_target, "w") as target_writer:
        target_writer.write(">" + hit.id + "\n" + str(hit.seq.translate()))

    temp_query = output + query + "_" + str(local_index) + "_query.faa"
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query_fasta[query].id + "\n" + str(query_fasta[query].seq))

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
                        length = int(line.split(" ")[2].strip())
                    if(line.startswith("# Identity")):
                        identity = float(line.split("(")[1][:-3].strip())
                    if(line.startswith("# Similarity")):
                        similarity = float(line.split("(")[1][:-3].strip())
                        break

    os.remove(temp_target)
    os.remove(temp_output)
    os.remove(temp_query)
    st_result = []
    hssp_identity = calculateHSSPIdentity(length, distance)
    hssp_similarity = calculateHSSPSimilarity(length, distance)
    if(similarity > identity and (identity >= hssp_identity or (use_sim > 0 and similarity >= hssp_similarity))):
        st_result = [id, subject, query, str(identity), str(similarity), str(length), str(hssp_identity), str(hssp_similarity)]
    print("\tTarget: " + hit.id + "\tfinished")
    return "\t".join(list(filter(None, st_result)))


def calculateHSSPIdentity(length, distance):
    if(length < 12):
        return 100
    elif(length < 418):
        return distance + 480 * length**(-0.32 * (1 + math.exp(-length / 1000)))
    else:
        return distance + 19.5

def calculateHSSPSimilarity(length, distance):
    if(length < 9):
        return 100
    elif(length < 742):
        return distance + 420 * length**(-0.335 * (1 + math.exp(-length / 2000)))
    else:
        return distance + 9.9

def mergeResults(config):
    for subject in config["subjects"]:
        for query in config["queries"]:
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
                ex_fasta = "alignment/exonerate/" + subject + "/" + query + ".cds"
                ex_stretcher = "evaluation/exonerate/" + subject + "/" + query + ".stretcher"
                ex_id_list = []
                sp_gff = "alignment/spaln/" + subject + "/" + query + ".gff"
                sp_fasta = "alignment/spaln/" + subject + "/" + query + ".cds"
                sp_stretcher = "evaluation/spaln/" + subject + "/" + query + ".stretcher"
                sp_id_list = []
                hit_id = 0
                if(not (os.path.exists(output + query + ".gff") or os.path.exists(output + query + ".cds"))):
                    if(os.path.exists(ex_gff) and os.path.exists(ex_stretcher) and os.path.exists(sp_gff) and os.path.exists(sp_stretcher)):
                        sp_gff_df = pd.read_csv(sp_gff, sep="\t", header=0)
                        sp_eval_df = pd.read_csv(sp_stretcher, sep="\t", header=0)
                        ex_gff_df = pd.read_csv(ex_gff, sep="\t", header=0)
                        ex_eval_df = pd.read_csv(ex_stretcher, sep="\t", header=0)
                        for sp_index,sp_row in sp_eval_df.iterrows():
                            if(not sp_row.empty):
                                ex_gene_rows = ex_gff_df.loc[(ex_gff_df["type"] == "gene") & (ex_gff_df["#subject"] == sp_row["target"]) & (ex_gff_df["query"] == sp_row["query"])]
                                ex_eval_rows = ex_eval_df.loc[(ex_eval_df["target"] == sp_row["target"]) & (ex_eval_df["query"] == sp_row["query"])]
                                ex_sequences = SeqIO.index(ex_fasta, "fasta")
                                sp_gene_rows = sp_gff_df.loc[(sp_gff_df["type"] == "gene") & (sp_gff_df["#subject"] == sp_row["target"]) & (sp_gff_df["query"] == sp_row["query"])]
                                sp_sequences = SeqIO.index(sp_fasta, "fasta")
                                sp_id_row = sp_gene_rows.loc[sp_gene_rows["id"] == str(sp_row["#id"])].iloc[0]
                                sp_header = sp_row["query"] + "::" + sp_row["target"] + "::" + str(sp_row["#id"])
                                for ex_index,ex_row in ex_eval_rows.iterrows():
                                    if(not ex_row.empty):
                                        if(not str(ex_row["#id"]) in ex_id_list):
                                            ex_id_row = ex_gene_rows.loc[ex_gene_rows["id"] == str(ex_row["#id"])].iloc[0]
                                            ex_header = ex_row["query"] + "::" + ex_row["target"] + "::" + str(ex_row["#id"])
                                            if(sp_row["target"] == ex_row["target"] and sp_row["query"] == ex_row["query"] and sp_id_row["orientation"] == ex_id_row["orientation"]):
                                                if(calculateOverlapPercentage(sp_id_row["start"], sp_id_row["end"], ex_id_row["start"], ex_id_row["end"]) >= config["overlap_percentage"]):
                                                    sp_id_list.append(str(sp_row["#id"]))
                                                    ex_id_list.append(str(ex_row["#id"]))
                                                    if((str(sp_sequences[sp_header].seq).startswith("ATG") and str(ex_sequences[ex_header].seq).startswith("ATG"))
                                                       or not (str(sp_sequences[sp_header].seq).startswith("ATG") or str(ex_sequences[ex_header].seq).startswith("ATG"))):
                                                       if(sp_id_row["score"] >= ex_id_row["score"]):
                                                           id_rows = sp_gff_df.loc[sp_gff_df["id"] == str(sp_row["#id"])]
                                                           id_rows["id"] = hit_id
                                                           header = "##" +  sp_row["target"] + "\t" + sp_row["query"] + "\t" + str(hit_id)
                                                           gff_list.append(header + "\n" + id_rows.to_string(header=False, index=False))
                                                           cds_list.append(">" + sp_row["query"] + "::" + sp_row["target"] + "::" + str(hit_id) + "\n" + str(sp_sequences[sp_header].seq))
                                                           hit_id += 1
                                                       else:
                                                           id_rows = ex_gff_df.loc[ex_gff_df["id"] == str(ex_row["#id"])]
                                                           id_rows["id"] = hit_id
                                                           header = "##" +  ex_row["target"] + "\t" + ex_row["query"] + "\t" + str(hit_id)
                                                           gff_list.append(header + "\n" + id_rows.to_string(header=False, index=False))
                                                           cds_list.append(">" + ex_row["query"] + "::" + ex_row["target"] + "::" + str(hit_id) + "\n" + str(ex_sequences[ex_header].seq))
                                                           hit_id += 1
                                                    elif(str(sp_sequences[sp_header].seq).startswith("ATG")):
                                                        id_rows = sp_gff_df.loc[sp_gff_df["id"] == str(sp_row["#id"])]
                                                        id_rows["id"] = hit_id
                                                        header = "##" +  sp_row["target"] + "\t" + sp_row["query"] + "\t" + str(hit_id)
                                                        gff_list.append(header + "\n" + id_rows.to_string(header=False, index=False))
                                                        cds_list.append(">" + sp_row["query"] + "::" + sp_row["target"] + "::" + str(hit_id) + "\n" + str(sp_sequences[sp_header].seq))
                                                        hit_id += 1
                                                    else:
                                                        id_rows = ex_gff_df.loc[ex_gff_df["id"] == str(ex_row["#id"])]
                                                        id_rows["id"] = hit_id
                                                        header = "##" +  ex_row["target"] + "\t" + ex_row["query"] + "\t" + str(hit_id)
                                                        gff_list.append(header + "\n" + id_rows.to_string(header=False, index=False))
                                                        cds_list.append(">" + ex_row["query"] + "::" + ex_row["target"] + "::" + str(hit_id) + "\n" + str(ex_sequences[ex_header].seq))
                                                        hit_id += 1
                                                    break

                        sp_ids = sp_eval_df.loc[sp_eval_df["#id"]]
                        for sp_index,sp_row in sp_ids.iterrows():
                            if(not (sp_row.empty or str(sp_row["#id"]) in sp_id_list)):
                                id_rows = sp_gff_df.loc[sp_gff_df["id"] == str(sp_row["#id"])]
                                id_rows["id"] = hit_id
                                if(not id_rows.empty and str(sp_row["target"]) != "nan" and str(sp_row["query"]) != "nan"):
                                    header = "##" +  sp_row["target"] + "\t" + sp_row["query"] + "\t" + str(hit_id)
                                    gff_list.append(header + "\n" + id_rows.to_string(header=False, index=False))
                                    cds_list.append(">" + sp_row["query"] + "::" + sp_row["target"] + "::" + str(hit_id) + "\n" + str(sp_sequences[sp_header].seq))
                                    hit_id += 1

                        ex_ids = ex_eval_df.loc[ex_eval_df["#id"]]
                        for ex_index,ex_row in ex_ids.iterrows():
                            if(not (ex_row.empty or str(ex_row["#id"]) in ex_id_list)):
                                id_rows = ex_gff_df.loc[ex_gff_df["id"] == str(ex_row["#id"])]
                                id_rows["id"] = hit_id
                                if(not id_rows.empty and str(ex_row["target"]) != "nan" and str(ex_row["query"]) != "nan"):
                                    header = "##" +  ex_row["target"] + "\t" + ex_row["query"] + "\t" + str(hit_id)
                                    gff_list.append(header + "\n" + id_rows.to_string(header=False, index=False))
                                    cds_list.append(">" + ex_row["query"] + "::" + ex_row["target"] + "::" + str(hit_id) + "\n" + str(ex_sequences[ex_header].seq))
                                    hit_id += 1

                    with open(output + query + ".gff", "w") as gff_writer:
                        gff_writer.write("\n".join(gff_list))
                    with open(output + query + ".cds", "w") as cds_writer:
                        cds_writer.write("\n".join(cds_list))

                    del gff_list[:]
                    del cds_list[:]
                    del ex_id_list[:]
                    del sp_id_list[:]
                    print("Merging: finished")
                else:
                    print("Result files found, skipping")
            except:
                print("Error while merging results, see log file '" + log_path + query + ".log'")
                with open(log_path + query + ".log", "w") as log_writer:
                    log_writer.write(str(traceback.format_exc()))

def calculateOverlapPercentage(start1, end1, start2, end2):
    if(end1 <= start2 or end2 < start1):
        return 0
    else:
        return ((min(int(end1), int(end2)) - max(int(start1), int(start2))) / min(int(end1)-int(start1), int(end2)-int(start2))) * 100



config = readConfigFile("config.yaml")
buildBLASTDatabase(config)
runBLASTCommand(config)
manager = mp.Manager()
index = manager.Value("i", -1)
lock = manager.Lock()
runExonerateCommand(config)
index.value = -1
runSpalnCommand(config)
index.value = -1
evaluateAlignment(config, "exonerate")
index.value = -1
evaluateAlignment(config, "spaln")
mergeResults(config)
