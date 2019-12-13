from Bio.Seq import Seq


def overlap(start1, end1, orientation1, start2, end2, orientation2):
    if(orientation1 != orientation2):
        return 0

    if(start1 <= start2 and end1 >= end2):
        return 100

    if(start1 >= start2 and end1 <= end2):
        return 100

    return ((min(end1, end2) - max(start1, start2)) / min(end1-start1, end2-start2)) * 100


rule Merge_Results:
    input:
        "scores/exonerate/{subject}/{query}.sc_gff",
        "scores/spaln/{subject}/{query}.sc_gff"
    output:
        "merged/{subject}/{query}.seq_gff"
    params:
        config["overlap_percentage"]
    log:
        "log/merged/{subject}/{query}.log"
    run:
        try:
            ex_gff = {}
            if(os.stat(input[0]).st_size != 0):
                with open(input[0], "r") as ex_reader:
                    content = ex_reader.readlines()
                    infos = {}
                    contig = None
                    query = None
                    for line in content:
                        if(line):
                            if(line.startswith("#")):
                                infos = {}
                                splitted_line = line.split("\t")
                                contig = splitted_line[0][1:].strip()
                                query = splitted_line[1].strip()
                                infos["start"] = int(splitted_line[2].strip())
                                infos["end"] = int(splitted_line[3].strip())
                                infos["pos"] = []
                                infos["orientation"] = splitted_line[4].strip()
                                infos["identity"] = float(splitted_line[5].strip())
                                infos["similarity"] = float(splitted_line[6].strip())
                            elif(line.startswith(">nuc")):
                                infos["nuc"] = line.split(">nuc:")[1].strip()
                            elif(line.startswith(">cds")):
                                infos["cds"] = line.split(">cds:")[1].strip()
                            elif(line.startswith(">pep")):
                                infos["pep"] = line.split(">pep:")[1].strip()
                                infos["stop"] = [i for i, ltr in enumerate(infos["pep"]) if ltr == "*"]
                                if(not contig in ex_gff):
                                    ex_gff[contig] = {}

                                if(not query in ex_gff[contig]):
                                    ex_gff[contig][query] = [infos]
                                else:
                                    no_overlap = True
                                    index_remove_list = []
                                    for hit in ex_gff[contig][query]:
                                        if(overlap(start, end, infos["orientation"], hit["start"], hit["end"], hit["orientation"]) >= params[0]):
                                            no_overlap = False
                                            if((infos["pep"].startswith("M") and hit["pep"].startswith("M")) or
                                               not (infos["pep"].startswith("M") and not hit["pep"].startswith("M"))):
                                                if(infos["identity"] > hit["identity"] or
                                                   (infos["identity"] == hit["identity"] and infos["similarity"] > hit["similarity"])):
                                                    index_remove_list.append(ex_gff[contig][query].index(hit))
                                            elif(infos["pep"].startswith("M")):
                                                index_remove_list.append(ex_gff[contig][query].index(hit))

                                    if(no_overlap):
                                        ex_gff[contig][query].append(infos)
                                    elif(len(index_remove_list)):
                                        for index in index_remove_list:
                                            del ex_gff[contig][query][index]

                                        ex_gff[contig][query].append(infos)
                            else:
                                splitted_line = line.split("\t")
                                start = int(splitted_line[2].strip())
                                end = int(splitted_line[3].strip())
                                infos["pos"].append(str(start) + ";" + str(end))

            sp_gff = {}
            if(os.stat(input[1]).st_size != 0):
                with open(input[1], "r") as ex_reader:
                    content = ex_reader.readlines()
                    infos = {}
                    contig = None
                    query = None
                    for line in content:
                        if(line):
                            if(line.startswith("#")):
                                infos = {}
                                splitted_line = line.split("\t")
                                contig = splitted_line[0][1:].strip()
                                query = splitted_line[1].strip()
                                infos["start"] = int(splitted_line[2].strip())
                                infos["end"] = int(splitted_line[3].strip())
                                infos["pos"] = []
                                infos["orientation"] = splitted_line[4].strip()
                                infos["identity"] = float(splitted_line[5].strip())
                                infos["similarity"] = float(splitted_line[6].strip())
                            elif(line.startswith(">nuc")):
                                infos["nuc"] = line.split(">nuc:")[1].strip()
                            elif(line.startswith(">cds")):
                                infos["cds"] = line.split(">cds:")[1].strip()
                            elif(line.startswith(">pep")):
                                infos["pep"] = line.split(">pep:")[1].strip()
                                infos["stop"] = [i for i, ltr in enumerate(infos["pep"]) if ltr == "*"]
                                if(not contig in sp_gff):
                                    sp_gff[contig] = {}

                                if(not query in sp_gff[contig]):
                                    sp_gff[contig][query] = [infos]
                                else:
                                    no_overlap = True
                                    index_remove_list = []
                                    for hit in sp_gff[contig][query]:
                                        if(overlap(start, end, infos["orientation"], hit["start"], hit["end"], hit["orientation"]) >= params[0]):
                                            no_overlap = False
                                            if((infos["pep"].startswith("M") and hit["pep"].startswith("M")) or
                                               not (infos["pep"].startswith("M") or hit["pep"].startswith("M"))):
                                                if(infos["identity"] > hit["identity"] or
                                                   (infos["identity"] == hit["identity"] and infos["similarity"] > hit["similarity"])):
                                                    index_remove_list.append(sp_gff[contig][query].index(hit))
                                            elif(infos["pep"].startswith("M")):
                                                index_remove_list.append(sp_gff[contig][query].index(hit))

                                    if(no_overlap):
                                        sp_gff[contig][query].append(infos)
                                    elif(len(index_remove_list)):
                                        for index in index_remove_list:
                                            del sp_gff[contig][query][index]

                                        sp_gff[contig][query].append(infos)
                            else:
                                splitted_line = line.split("\t")
                                start = int(splitted_line[2].strip())
                                end = int(splitted_line[3].strip())
                                infos["pos"].append(str(start) + ";" + str(end))

            joined_results = []
            inner_joined_results = []
            for contig in ex_gff:
                if(contig in sp_gff):
                    for query in ex_gff[contig]:
                        if(query in sp_gff[contig]):
                            for ex_info in ex_gff[contig][query]:
                                no_overlap = True
                                for sp_info in sp_gff[contig][query]:
                                    if(overlap(ex_info["start"], ex_info["end"], ex_info["orientation"], sp_info["start"], sp_info["end"], sp_info["orientation"]) >= params[0]):
                                        no_overlap = False
                                        if((ex_info["pep"].startswith("M") and sp_info["pep"].startswith("M")) or
                                           not (ex_info["pep"].startswith("M") or sp_info["pep"].startswith("M"))):
                                            if(ex_info["identity"] > sp_info["identity"] or
                                               (ex_info["identity"] == sp_info["identity"] and ex_info["similarity"] > sp_info["similarity"])):
                                                if(not ex_info in inner_joined_results):
                                                    inner_joined_results.append(ex_info)
                                        elif(ex_info["pep"].startswith("M")):
                                            if(not ex_info in inner_joined_results):
                                                inner_joined_results.append(ex_info)

                                if(no_overlap):
                                    if(not ex_info in inner_joined_results):
                                        inner_joined_results.append(ex_info)

                            for sp_info in sp_gff[contig][query]:
                                no_overlap = True
                                for ex_info in ex_gff[contig][query]:
                                    if(overlap(sp_info["start"], sp_info["end"], sp_info["orientation"], ex_info["start"], ex_info["end"], ex_info["orientation"]) >= params[0]):
                                        no_overlap = False
                                        if(not ex_info in inner_joined_results):
                                            if((sp_info["pep"].startswith("M") and ex_info["pep"].startswith("M")) or
                                               not (sp_info["pep"].startswith("M") or ex_info["pep"].startswith("M"))):
                                                if(sp_info["identity"] > ex_info["identity"] or
                                                   (sp_info["identity"] == ex_info["identity"] and sp_info["similarity"] > ex_info["similarity"])):
                                                    if(not sp_info in inner_joined_results):
                                                        inner_joined_results.append(sp_info)
                                            elif(sp_info["pep"].startswith("M")):
                                                if(not sp_info in inner_joined_results):
                                                    inner_joined_results.append(sp_info)

                                if(no_overlap):
                                    if(not sp_info in inner_joined_results):
                                        inner_joined_results.append(sp_info)

                            for info in inner_joined_results:
                                joined_results.append(contig + "\t" + query + "\t" + str(info["start"]) + "\t" + str(info["end"]) +
                                                      "\t" + info["orientation"] + "\t" + "|".join(info["pos"]) + "\t" + str(info["identity"]) +
                                                      "\t" + str(info["similarity"]) + "\t" + info["nuc"] + "\t" + info["cds"] +
                                                      "\t" + info["pep"] + "\t" + "|".join(info["stop"]))

                            del inner_joined_results[:]
                        else:
                            for info in ex_gff[contig][query]:
                                joined_results.append(contig + "\t" + query + "\t" + str(info["start"]) + "\t" + str(info["end"]) + "\t" +
                                                      info["orientation"] + "\t" + "|".join(info["pos"]) + "\t" + str(info["identity"]) +
                                                      "\t" + str(info["similarity"]) + "\t" + info["nuc"] + "\t" + info["cds"] + "\t" + info["pep"]
                                                       + "\t" + "|".join(info["stop"])
                    for query in sp_gff[contig]:
                        if(not query in ex_gff[contig]):
                            for info in sp_gff[contig][query]:
                                joined_results.append(contig + "\t" + query + "\t" + str(info["start"]) + "\t" + str(info["end"]) + "\t" +
                                                      info["orientation"] + "\t" + "|".join(info["pos"]) + "\t" + str(info["identity"]) +
                                                      "\t" + str(info["similarity"]) + "\t" + info["nuc"] + "\t" + info["cds"] + "\t" + info["pep"]
                                                       + "\t" + "|".join(info["stop"])
                else:
                    for query in ex_gff[contig]:
                        for info in ex_gff[contig][query]:
                            joined_results.append(contig + "\t" + query + "\t" + str(info["start"]) + "\t" + str(info["end"]) + "\t" +
                                                  info["orientation"] + "\t" + "|".join(info["pos"]) + "\t" + str(info["identity"]) +
                                                  "\t" + str(info["similarity"]) + "\t" + info["nuc"] + "\t" + info["cds"] + "\t" + info["pep"]
                                                   + "\t" + "|".join(info["stop"])

            for contig in sp_gff:
                if(not contig in ex_gff):
                    for query in sp_gff[contig]:
                        for info in sp_gff[contig][query]:
                            joined_results.append(contig + "\t" + query + "\t" + str(info["start"]) + "\t" + str(info["end"]) + "\t" +
                                                  info["orientation"]+ "\t" + "|".join(info["pos"]) + "\t" + str(info["identity"]) +
                                                  "\t" + str(info["similarity"]) + "\t" + info["nuc"] + "\t" + info["cds"] + "\t" + info["pep"])

            with open(output[0], "w") as merged_writer:
                merged_writer.write("\n".join(joined_results))

            ex_gff.clear()
            sp_gff.clear()
            del joined_results[:]
            if(os.path.exists(output[0])):
                os.system("touch " + output[0])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()) + "\n" + str(ex))
