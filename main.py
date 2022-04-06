"""
Comp Bio Assignment 4
Written by Ray Loerke
BioGRID MITAB File Processor
"""


def main():
    process_biogrid_file("BIOGRID-ORGANISM-Homo_sapiens-4.4.207.mitab.txt", "H_Sapiens_BioGRID_network.txt")


def get_physical_interaction_codes(code_list):
    """
    This function finds all MI codes descending from and including MI: 0407 (direct interaction)

    :param code_list: This is a set of strings containing all MI codes and their is_a designations
    :return: Returns a set of MI codes that belong to the 0407 root
    """
    return mio_helper(code_list, 407, set())


def get_experimental_codes(code_list):
    """
    This function finds all MI codes descending from and including MI: 0045 (experimental detection methods)

    :param code_list: This is a set of strings containing all MI codes and their is_a designations
    :return: Returns a set of MI codes that belong to the 0045 root
    """
    return mio_helper(code_list, 45, set())


def mio_helper(code_list, root_code, cat_list):
    """
    This function recursively searches a set of all MI codes to find all the codes descending from the given root code

    :param code_list: This is a set of strings containing all MI codes and their is_a designations
    :param root_code: This is the MI code that will serve as our root, we will find all MI codes descending from this
    :param cat_list: This is a set of descending MI codes where we will build our output,
    when this function is first called an empty set should be input for this parameter
    :return: When no change is made to the cat_list recursion ends and the cat_list set is returned
    """

    # These bools are used to track if a new MI code has been found and if recursion should continue
    continue_search = False
    new_term = True

    # If this is the first time this function has been called we seed cat_list with our root node
    if not cat_list:
        cat_list.add(root_code)

    # We loop through every MI code
    for code in code_list:

        # For every code we need to check if it is a descendant of a code already in out cat_list
        for term_id in cat_list:

            # Some formatting is done to pad the is_a string with zeroes
            # and then the MI code entry is searched for a matching is_a string
            if code.find("is_a: MI:" + ("0" * (4 - len(str(term_id)))) + str(term_id)) != -1:

                # If we find a matching is_a string we need to check if we already have this MI code in our cat_list
                for term_id_2 in cat_list:
                    if ("0" * (4 - len(str(term_id_2))) + str(term_id_2)) == code[7:11]:
                        new_term = False

                # If this is a new term then we add it to our cat_list set and indicate that recursion should continue
                if new_term:
                    cat_list.add(int(code[7:11]))
                    continue_search = True
                    break

        # The bool is reset for the next code
        new_term = True

    # If we found a new code, continue recursion, otherwise, return the final cat_list
    if continue_search:
        return mio_helper(code_list, root_code, cat_list)
    else:
        return cat_list


def format_helper(filename):
    """
    This function generates a set of all MI codes and their is_a designations from a .owl text file

    :param filename: The path of the text file containing the MI ontology
    :return: A set of strings of all MI codes and their is_a designations
    """

    # code_list will be the set that we return,
    # buffer will be used to build the strings before they are added to the set,
    # in_term tracks if we are in or between terms
    code_list = set()
    buffer = ""
    in_term = False

    # Open the ontology file and loop through it line by line
    mio = open(filename, "r")
    for line in mio:

        # If we are in a term and hit a newline character than the current term has ended,
        # add the string to the code_list, reset the buffer, and indicate that the term has ended
        if in_term and line == "\n":
            code_list.add(buffer)
            buffer = ""
            in_term = False

        # If we are in a term and find an id or is_a line add it to the buffer
        if in_term:
            if line.startswith("id:") or line.startswith("is_a:"):
                buffer += line

        # If we are not in a term and see the [Term] tag then indicate that we are moving into a new term
        if not in_term and line == "[Term]\n":
            in_term = True

    # Once we have gone through the whole file close it and return our set
    mio.close()
    return code_list


def process_biogrid_file(filename, path):
    """
    This function builds a list of protein-protein interactions that fit a set of parameters
    Each interaction added into the network must:
    -Consist of 2 human (taxid: 9606) proteins
    -Be a physical interaction based on its MI code
    -Be experimentally detected based on its MI code
    -Not be a duplicate of another interaction
    -Not be a self-loop

    :param filename: The path to the text file containing the BioGRID Mitab text file
    :param path: The path to the network file to be generated by this function
    :return: Nothing is returned by this function,
    it generates a network text file containing the desired protein-protein interactions
    This file starts with a series of comment lines starting with '#' describing the statistics of network generation
    It then lists the interactions on separate lines, listed by their Entrez Gene names, separated by a tab
    The list is alphabetically ordered inter and intra interaction
    """

    # Statistic tracking variables are created and initialized
    interactions_read = 0
    interactions_kept = 0
    interactions_in_network = 0
    discarded_inter_species = 0
    discarded_physical = 0
    discarded_experimental = 0
    discarded_duplicate = 0
    discarded_self_loop = 0

    # Interactions will be the set that holds the valid interactions, proteins will track each unique protein
    interactions = set()
    proteins = set()

    # These bools will track if the interaction is valid and if we have moved past the BioGRID's header
    valid = True
    header = True

    # Use a helper function to generate a set of all MI codes and their is_a designations
    code_list = format_helper("mi.owl")

    # Pass this set to functions to generate new sets with only the MI codes for physical interaction
    # and experimentally detected interactions
    phys_set = get_physical_interaction_codes(code_list)
    exp_set = get_experimental_codes(code_list)

    # Open the BioGRID file and loop through it line by line
    bio = open(filename, "r")
    for line in bio:

        # Skip the line if it is the header line
        if header:
            header = False
        else:

            # Split the line into a list delineated by the tab character,
            # split the Alt ID sections further to isolate the Entrez Gene names
            interactions_read += 1
            line_list = line.split('\t')
            id_a_list = line_list[2].split('|')
            id_b_list = line_list[3].split('|')

            # If the interaction is not between human proteins reject it
            if line_list[9] != "taxid:9606" or line_list[10] != "taxid:9606":
                discarded_inter_species += 1
                valid = False

            # If the interactions MI code is not in the set of physical interactions reject it
            elif not int(line_list[11][11:15]) in phys_set:
                discarded_physical += 1
                valid = False

            # If the interactions MI code is not in set of experimentally detected interactions reject it
            elif not int(line_list[6][11:15]) in exp_set:
                discarded_experimental += 1
                valid = False

            # If the interaction is a self-loop discard it
            elif line_list[0] == line_list[1]:
                discarded_self_loop += 1
                valid = False

            # If the interaction is a duplicate reject it
            elif (id_a_list[1][22:] + "\t" + id_b_list[1][22:]) in interactions:
                discarded_duplicate += 1
                valid = False

            # Check for duplicate interactions ordered the opposite way
            elif (id_b_list[1][22:] + "\t" + id_a_list[1][22:]) in interactions:
                discarded_duplicate += 1
                valid = False

            # If the interaction has passed these tests it is a valid interaction
            if valid:
                interactions_kept += 1

                # Add the interaction to the set with the interactors ordered alphabetically
                if id_a_list[1][22:] < id_b_list[1][22:]:
                    interactions.add(id_a_list[1][22:] + "\t" + id_b_list[1][22:])
                else:
                    interactions.add(id_b_list[1][22:] + "\t" + id_a_list[1][22:])

                # If the first protein in the interaction is new add it to the protein list
                if not id_a_list[1][22:] in proteins:
                    interactions_in_network += 1
                    proteins.add(id_a_list[1][22:])

                # If the second protein in the interaction is new add it to the protein list
                if not id_b_list[1][22:] in proteins:
                    interactions_in_network += 1
                    proteins.add(id_b_list[1][22:])

            # Reset the bool
            valid = True

    # Once all interactions have been read close the file
    bio.close()

    # Start the output string
    out_string = ""

    # Add the statistics headers to the output string
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

    # Add the interactions to the output string
    for item in sorted(interactions):
        out_string += item
        out_string += "\n"

    # Write the output string to a file using the given path
    out_file = open(path, "w")
    out_file.write(out_string)
    out_file.close()

    # Notify the user that the network file has been created
    print("Network File Created!")


main()
