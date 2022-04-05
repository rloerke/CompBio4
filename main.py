"""
Comp Bio Assignment 4
Written by Ray Loerke
BioGRID MITAB File Processor
"""


def main():
    """
    code_list = format_helper("mi.owl")

    phys_set = get_physical_interaction_codes(code_list)
    print(phys_set)
    print(len(phys_set))
    exp_set = get_experimental_codes(code_list)
    print(exp_set)
    print(len(exp_set))
    """

    process_biogrid_file("BIOGRID-ORGANISM-Homo_sapiens-4.4.207.mitab.txt", "H_Sapiens_BioGRID_network.txt")


def get_physical_interaction_codes(code_list):
    return mio_helper(code_list, 407, set())


def get_experimental_codes(code_list):
    return mio_helper(code_list, 45, set())


def mio_helper(code_list, root_code, cat_list):
    continue_search = False
    new_term = True

    if not cat_list:
        cat_list = set()
        cat_list.add(root_code)

    for code in code_list:
        for term_id in cat_list:
            if code.find("is_a: MI:" + ("0" * (4 - len(str(term_id)))) + str(term_id)) != -1:
                for term_id_2 in cat_list:
                    if ("0" * (4 - len(str(term_id_2))) + str(term_id_2)) == code[7:11]:
                        new_term = False
                if new_term:
                    cat_list.add(int(code[7:11]))
                    continue_search = True
                    break
        new_term = True

    if continue_search:
        return mio_helper(code_list, root_code, cat_list)
    else:
        return cat_list


def format_helper(filename):
    code_list = set()
    buffer = ""
    in_term = False

    mio = open(filename, "r")
    for line in mio:

        if in_term and line == "\n":
            code_list.add(buffer)
            buffer = ""
            in_term = False

        if in_term:
            if line.startswith("id:") or line.startswith("is_a:"):
                buffer += line

        if not in_term and line == "[Term]\n":
            in_term = True

    mio.close()
    return code_list


def process_biogrid_file(filename, path):
    interactions_read = 0
    interactions_kept = 0
    interactions_in_network = 0
    discarded_inter_species = 0
    discarded_physical = 0
    discarded_experimental = 0
    discarded_duplicate = 0
    discarded_self_loop = 0

    interactions = set()
    proteins = set()
    valid = True
    header = True

    code_list = format_helper("mi.owl")
    phys_set = get_physical_interaction_codes(code_list)
    exp_set = get_experimental_codes(code_list)

    bio = open(filename, "r")

    for line in bio:
        if header:
            header = False
        else:
            interactions_read += 1
            line_list = line.split('\t')
            id_a_list = line_list[2].split('|')
            id_b_list = line_list[3].split('|')

            if line_list[9] != "taxid:9606" or line_list[10] != "taxid:9606":
                discarded_inter_species += 1
                valid = False

            elif not int(line_list[11][11:15]) in phys_set:
                discarded_physical += 1
                valid = False

            elif not int(line_list[6][11:15]) in exp_set:
                discarded_experimental += 1
                valid = False

            elif line_list[0] == line_list[1]:
                discarded_self_loop += 1
                valid = False

            elif (id_a_list[1][22:] + "\t" + id_b_list[1][22:]) in interactions:
                discarded_duplicate += 1
                valid = False

            elif (id_b_list[1][22:] + "\t" + id_a_list[1][22:]) in interactions:
                discarded_duplicate += 1
                valid = False

            if valid:
                interactions_kept += 1

                if id_a_list[1][22:] < id_b_list[1][22:]:
                    interactions.add(id_a_list[1][22:] + "\t" + id_b_list[1][22:])
                else:
                    interactions.add(id_b_list[1][22:] + "\t" + id_a_list[1][22:])

                if not id_a_list[1][22:] in proteins:
                    interactions_in_network += 1
                    proteins.add(id_a_list[1][22:])

                if not id_b_list[1][22:] in proteins:
                    interactions_in_network += 1
                    proteins.add(id_b_list[1][22:])

            valid = True

    bio.close()

    out_string = ""

    out_string += "#Interactions Read:\t" + str(interactions_read) + "\n"
    out_string += "#Interactions Kept In Network:\t" + str(interactions_kept) + "\t" + \
                  str(round(interactions_kept / interactions_read * 100, 2)) + "%" + "\n"
    out_string += "#Total Number of Proteins:\t" + str(interactions_in_network) + "\n"
    out_string += "#Discarded Inter-Species:\t" + str(discarded_inter_species) + "\t" + \
                  str(round(discarded_inter_species / interactions_read * 100, 2)) + "%" + "\n"
    out_string += "#Discarded Not Physical:\t" + str(discarded_physical) + "\t" + \
                  str(round(discarded_physical / interactions_read * 100, 2)) + "%" + "\n"
    out_string += "#Discarded Not Experimental:\t" + str(discarded_experimental) + "\t" + \
                  str(round(discarded_experimental / interactions_read * 100, 2)) + "%" + "\n"
    out_string += "#Discarded Self-Loop:\t" + str(discarded_self_loop) + "\t" + \
                  str(round(discarded_self_loop / interactions_read * 100, 2)) + "%" + "\n"
    out_string += "#Discarded Duplicate:\t" + str(discarded_duplicate) + "\t" + \
                  str(round(discarded_duplicate / interactions_read * 100, 2)) + "%" + "\n"

    for item in sorted(interactions):
        out_string += item
        out_string += "\n"

    out_file = open(path, "w")
    out_file.write(out_string)
    out_file.close()
    
    print("Network File Created!")


main()
