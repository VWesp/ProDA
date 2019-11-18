import os

for subject in snakemake.input:
    if(subject.endswith(".fasta")):
        subject_name = subject.split("/")[-1].split(".fasta")[0]
        sequence = []
        subject_content = {}
        with open(subject, "r") as subject_reader:
            content = subject_reader.readlines()
            current_header = ""
            for line in content:
                if(line.startswith(">")):
                    if(not subject_content):
                        subject_content[line[1:].strip()] = None
                    else:
                        subject_content[current_header] = "".join(sequence)
                        del sequence[:]

                    current_header = line[1:].strip()
                else:
                    sequence.append(line.strip())

            subject_content[current_header] = "".join(sequence)
            del sequence[:]

        '''query_content = {}
        with open(snakemake.input[1], "r") as query_reader:
            content = query_reader.readlines()
            current_header = ""
            for line in content:
                if(line.startswith(">")):

                    if(not query_content):
                        query_content[line[1:].strip()] = None
                    else:
                        query_content[current_header] = "".join(sequence)
                        del sequence[:]

                    current_header = line[1:].strip()
                else:
                    sequence.append(line.strip())

            query_content[current_header] = "".join(sequence)
            del sequence[:]'''

        for hit in snakemake.input:
            if(hit.endswith(".hit")):
                hit_name = hit.split("/")[-2]
                if(hit_name == subject_name):
                    query_subject_matcher = {}
                    subject_start = None
                    subject_end = None
                    with open(hit, "r") as blast_reader:
                        content = blast_reader.readlines()
                        for line in content:
                            line_splitted = line.split("\t")
                            query_id = line_splitted[0]
                            subject_id = line_splitted[1]
                            subject_start = int(line_splitted[2])
                            subject_end = int(line_splitted[3])
                            if(not query_id in query_subject_matcher):
                                query_subject_matcher[query_id] = {}

                            if(not subject_id in query_subject_matcher[query_id]):
                                query_subject_matcher[query_id][subject_id] = []

                            if(subject_start > subject_end):
                                temp_loc = subject_start
                                subject_start = subject_end
                                subject_end = temp_loc

                            if(subject_start - snakemake.params[0] < 0):
                                subject_start = 0
                            else:
                                subject_start -= snakemake.params[0]

                            if(subject_end + snakemake.params[1] >= len(subject_content[subject_id])):
                                subject_end = len(subject_content[subject_id])-1
                            else:
                                subject_end += snakemake.params[1]

                            sequence = subject_content[subject_id][subject_start:subject_end]
                            query_subject_matcher[query_id][subject_id].append(sequence)

                    for query_id in query_subject_matcher:
                        for subject_id, sequences in query_subject_matcher[query_id].items():
                            query_indexer = 0
                            for sequence in sequences:
                                path = "matches/" + subject_name + "/" + query_id
                                if(not os.path.exists(path)):
                                    os.makedirs(path)

                                file_name = subject_id + "_start:" + str(subject_start) + "_end:" + str(subject_end) + "__" + str(query_indexer)
                                with open(path + "/" + subject_id + ".fasta", "a") as match_writer:
                                    match_writer.write(">" + file_name + "\n" + sequence + "\n")

                                query_indexer += 1

os.system("touch temp/scripts.txt")
