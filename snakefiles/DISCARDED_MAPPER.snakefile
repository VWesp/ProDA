import os
import shutil
from Bio.Seq import Seq
import traceback
import multiprocessing as mp


def mapResultsMultiprocessing(line, subjects, subject_name):
    map_results = []
    ex_results = []
    pos_list = []
    nuc_list = []
    stripped_line = line.strip()
    if(stripped_line):
        subject = stripped_line.split("\t")[0]
        query = stripped_line.split("\t")[1]
        idenity = line.split("\t")[3].strip()
        similarity = line.split("\t")[4].strip()
        aa_sequence = stripped_line.split("\t")[5]
        query_sequence = stripped_line.split("\t")[6]
        if(subject == subject_name):
            hit_output =  "results/discarded/mapped/" + subject + "/" + query + ".faa"
            with open(hit_output, "w") as hit_writer:
                hit_writer.write(">" + query + "\n" + aa_sequence)

            log_output = "log/results/discarded/mapped/" + subject + "/" + query + ".log"
            del ex_results[:]
            for fasta in SeqIO.parse(open(subjects), "fasta"):
                subject_output = "results/discarded/mapped/" + subject + "/" + query + "_subject.faa"
                with open(subject_output, "w") as hit_writer:
                    hit_writer.write(">" + fasta.id + "\n" + str(fasta.seq))

                map_output = "results/discarded/mapped/" + subject + "/" + query + ".ryo"
                os.system("(exonerate --model protein2genome --targettype dna --querytype protein "
                          "--ryo '>%ti::sstart=%tcb::send=%tce\n%tcs' --showalignment no --showvulgar no "
                          "--refine region --percent 100 "
                          "--query " + hit_output + " --target " + subject_output +
                          " > " + map_output + ") 2> " + log_output)

                seq_output = "results/discarded/mapped/" + subject + "/" + query + ".fna"
                os.system("(tail -n +4 " + map_output + " | head -n -1) > " + seq_output)
                with open(seq_output, "r") as seq_reader:
                    ex_results.append(seq_reader.read().strip())

                os.remove(subject_output)
                os.remove(map_output)
                os.remove(seq_output)

            os.remove(hit_output)
            for hits in ex_results:
                if(len(hits)):
                    for hit in hits.split(">"):
                        lines = hit.split("\n")
                        contig = lines[0].split("::")[0].strip()
                        sstart = lines[0].split("sstart=")[-1].split("::")[0].strip()
                        send = lines[0].split("send=")[-1].strip()
                        pos_list.append(contig + ":{" + sstart + "," + send + "}")
                        nuc_sequence = "".join(lines[1:]).strip()
                        if(not nuc_sequence in nuc_list):
                            nuc_list.append(nuc_sequence)

            map_results.append(subject + "\t" + query + "\t" + ";".join(pos_list) + "\t" + idenity + "\t" + similarity + "\t" +
                               ";".join(nuc_list) + "\t" + aa_sequence + "\t" + query_sequence)
            del pos_list[:]
            del nuc_list[:]

    return list(filter(None, map_results))



rule Map_Discarded_Results:
    input:
        "data/subjects/{subject}.fna",
        "results/discarded/proda_temp.tsv"
    output:
        temp("results/discarded/proda_{subject}.tsv")
    threads: config["threads"]
    run:
        try:
            if(os.stat(input[1]).st_size != 0):
                subject = input[0].split("/")[-1].split(".fna")[0]
                os.makedirs("results/discarded/mapped/" + subject)
                os.makedirs("log/results/discarded/mapped/" + subject)
                with open(input[1], "w") as discarded_reader:
                    lines = discarded_reader.read().split("\n")[1:]
                    pool = mp.Pool(processes=threads)
                    pool_map = partial(mapResultsMultiprocessing, subjects=input[0], subject_name=subject)
                    map_mp_results = pool.map_async(mapResultsMultiprocessing, lines)
                    pool.close()
                    pool.join()
                    map_joined_results = [r.strip() for result in map_mp_results.get() for r in result]

                with open(output[0], "w") as result_writer:
                    result_writer.write("Subject\tQuery\tPositions\tIdentity\tSimilarity\tNucleotide sequence\tAminoacid sequence\tQuery sequence\n"
                                        + "\n".join(map_joined_results))

                del map_joined_results[:]
                del map_mp_results.get()[:]
                shutil.rmtree("results/discarded/mapped/" + subject)
        except:
            print("\033[1;31;mError: " + str(traceback.format_exc()) + "\nSee log file: log/results/discarded/mapped.log")
            if(not os.path.exists("log/results/discarded/mapped")):
                os.makedirs("log/results/discarded/mapped")

            with open("log/results/discarded/mapped.log", "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))

        os.system("touch map_discarded_finished.txt")
