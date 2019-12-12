from Bio.Seq import Seq


configfile: "config.yaml"


def readGFF(gff, subjects):
    gff_list = []
    sequence = []
    last_orientation = None
    for result in gff:
        lines = result.split("\n")
        for line in lines:
            if(line and not line.startswith("#")):
                splitted_line = line.split("\t")
                contig = splitted_line[0].split("_cstart:")[0].strip()
                contig_start = int(splitted_line[0].split("_cstart:")[-1].split("_cend:")[0].strip())
                query = splitted_line[0].split("_query:")[-1].strip()
                type = splitted_line[2]
                start = int(splitted_line[3].strip()) + contig_start - 1
                end = int(splitted_line[4].strip()) + contig_start
                orientation = splitted_line[6].strip()
                if(type == "gene"):
                    if(len(sequence)):
                        seq = None
                        if(last_orientation == "-"):
                            seq = Seq("".join(sequence[::-1]))
                        else:
                            seq = Seq("".join(sequence))

                        gff_list.append(">cds:" + str(seq))
                        gff_list.append(">pep:" + str(seq.translate()))

                    last_orientation = orientation
                    gff_list.append("#" + contig + "\t" + query + "\t" + str(start) + "\t" + str(end) + "\t" + orientation)
                    if(orientation == "-"):
                        gff_list.append(">nuc:" + str(subjects[contig].seq[start:end].reverse_complement()).upper())
                    else:
                        gff_list.append(">nuc:" + str(subjects[contig].seq[start:end]).upper())

                    del sequence[:]
                elif(type == "cds"):
                    if(orientation == "-"):
                        sequence.append(str(subjects[contig].seq[start:end].reverse_complement()).upper())
                    else:
                        sequence.append(str(subjects[contig].seq[start:end]).upper())

                    gff_list.append(contig + "\t" + query + "\t" + str(start) + "\t" + str(end) + "\t" + orientation)

    if(len(sequence)):
        seq = None
        if(last_orientation == "-"):
            seq = Seq("".join(sequence[::-1]))
        else:
            seq = Seq("".join(sequence))

        gff_list.append(">cds:" + str(seq))
        gff_list.append(">pep:" + str(seq.translate()))

    return gff_list


rule finish:
    input:
        "finished.txt"

include: "snakefiles/BLAST.snakefile"

include: "snakefiles/QUERY_CONTIG_MATCHER.snakefile"

include: "snakefiles/EXONERATE.snakefile"

include: "snakefiles/EX_STRETCHER.snakefile"

include: "snakefiles/SPALN.snakefile"

include: "snakefiles/SP_STRETCHER.snakefile"

include: "snakefiles/JOIN_RESULTS.snakefile"

include: "snakefiles/THRESHOLD_FILTER.snakefile"

include: "snakefiles/VISUALIZE_RETAINED_RESULTS.snakefile"

include: "snakefiles/VISUALIZE_DISCARDED_RESULTS.snakefile"

rule finished:
    input:
        expand("temp/above/{subject}/{query}.txt", subject=config["subjects"], query=config["queries"]),
        expand("temp/below/{subject}/{query}.txt", subject=config["subjects"], query=config["queries"]),
    output:
        temp("finished.txt")
    shell:
        "touch {output} && rm -r temp"
