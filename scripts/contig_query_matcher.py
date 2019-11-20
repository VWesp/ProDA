import os
from Bio import SeqIO

try:
    subjects = [s for s in snakemake.input if s.endswith(".fna")]
    hits = [h for h in snakemake.input if h.endswith(".hit")]
    subjects_contigs = {}
    for subject in subjects:
        subject_name = subject.split("/")[-1].split(".fna")[0]
        subjects_contigs[subject_name] = {}
        fasta_sequences = SeqIO.parse(open(subject), "fasta")
        for fasta in fasta_sequences:
            subjects_contigs[subject_name][fasta.id] = str(fasta.seq)

    query_subject_matcher = {}
    for hit in hits:
        hit_name = hit.split("/")[-1].split(".hit")[0]
        subject_name = hit.split("/")[-2]
        if(not subject_name in query_subject_matcher):
            query_subject_matcher[subject_name] = {}

        with open(hit, "r") as blast_reader:
            content = blast_reader.readlines()
            for line in content:
                line_splitted = line.split("\t")
                contig_id = line_splitted[1]
                hit_start = int(line_splitted[2])
                hit_end = int(line_splitted[3])
                if(not hit_name in query_subject_matcher[subject_name]):
                    query_subject_matcher[subject_name][hit_name] = {}

                if(not contig_id in query_subject_matcher[subject_name][hit_name]):
                    query_subject_matcher[subject_name][hit_name][contig_id] = []

                if(hit_start > hit_end):
                    temp_loc = hit_start
                    hit_start = hit_end
                    hit_end = temp_loc

                if(hit_start - snakemake.params[0] < 0):
                    hit_start = 0
                else:
                    hit_start -= snakemake.params[0]

                if(hit_end + snakemake.params[1] >= len(subjects_contigs[subject_name][contig_id])):
                    hit_end = len(subjects_contigs[subject_name][contig_id])-1
                else:
                    hit_end += snakemake.params[1]

                sequence = subjects_contigs[subject_name][contig_id][hit_start:hit_end]
                query_subject_matcher[subject_name][hit_name][contig_id].append([sequence, hit_start, hit_end])

    subjects_contigs.clear()
    fasta = []
    for subject in query_subject_matcher:
        for query in query_subject_matcher[subject]:
            path = "matches/" + subject
            if(not os.path.exists(path)):
                os.makedirs(path)

            del fasta[:]
            for contig, sequences in query_subject_matcher[subject][query].items():
                for sequence in sequences:
                    header = contig + "_start:" + str(sequence[1]) + "_end:" + str(sequence[2])
                    fasta.append(">" + header + "\n" + sequence[0])

            with open(path + "/" + query + ".fna", "w") as match_writer:
                match_writer.write("\n".join(fasta))

    del fasta[:]
    query_subject_matcher.clear()
except Exception as ex:
    with open(snakemake.log[0], "w") as log_writer:
        log_writer.write(ex)
