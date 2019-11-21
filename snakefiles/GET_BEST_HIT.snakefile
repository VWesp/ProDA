from Bio import SeqIO, pairwise2

rule get_best_hit:
    input:
        "scores/{subject}/{query}.sc",
        "cd_hit/{subject}/{query}_merged.faa"
    output:
        "best_hit/{subject}/{query}.faa"
    log:
        "log/best_hit/{subject}/{query}.log"
    run:
        with open(input[0], "r") as score_reader:
            content = score_reader.readlines()
            scores = {}
            current_query = None
            index = None
            for line in content:
                if(line.strip()):
                    if(line.startswith("#")):
                        current_query = line[1:].strip()
                        scores[current_query] = ["", -1, -1]
                        index = 0
                    else:
                        contig = line.split("\t")[0].strip()
                        score = float(line.split("\t")[-1].strip())
                        if(scores[current_query][0] == contig):
                            index += 1

                        if(score > scores[current_query][1]):
                            scores[current_query] = [contig, score, index]

        best_hits = []
        fasta_sequences = SeqIO.parse(open(input[1]), "fasta")
        for fasta in fasta_sequences:
            query = fasta.id.split("::query=")[-1]
            index = 0
            if(scores[query][0] == fasta.id):
                if(scores[query][2] == index):
                    id = str(float(scores[query][1]) / len(str(fasta.seq)))
                    header = ">" + fasta.id + "::identity=" + id
                    best_hits.append(header + "\n" + str(fasta.seq))
                else:
                    index += 1

        with open(output[0], "w") as best_hit_writer:
            best_hit_writer.write("\n".join(best_hits))

        scores.clear()
        del best_hits[:]
