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

    process_biogrid_file("BIOGRID-ORGANISM-Homo_sapiens-4.4.207.mitab.txt", "network.txt")


def get_physical_interaction_codes(code_list):
    return mio_helper(code_list, 915, set())


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

            id_a_list = line_list[2].split('|')
            id_b_list = line_list[3].split('|')

            if (id_a_list[1][22:] + "\t" + id_b_list[1][22:]) in interactions:
                discarded_duplicate += 1
                valid = False

            if valid:
                interactions_kept += 1
                interactions.add(id_a_list[1][22:] + "\t" + id_b_list[1][22:])

            valid = True

    bio.close()

    print("\nInteractions Read: " + str(interactions_read))
    print("Interactions Kept In Network: " + str(interactions_kept))
    print("Total Interactors: " + str(interactions_in_network))
    print("Discarded Inter-Species: " + str(discarded_inter_species))
    print("Discarded Not Physical: " + str(discarded_physical))
    print("Discarded Not Experimental: " + str(discarded_experimental))
    print("Discarded Self-Loop: " + str(discarded_self_loop))
    print("Discarded Duplicate: " + str(discarded_duplicate))

    """
    out_string = ""

    for item in sorted(interactions):
        out_string += item
        out_string += "\n"

    out_file = open(path, "w")
    out_file.write(out_string)
    out_file.close()
    
    print("Network File Created!")
    """


main()
