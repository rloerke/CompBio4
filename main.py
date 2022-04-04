"""
Comp Bio Assignment 4
Written by Ray Loerke
BioGRID MITAB File Processor
"""


def main():
    print("Hello World!")
    get_physical_interaction_codes("mi.owl")


def get_physical_interaction_codes(filename):
    print("Phys Codes")
    mio_fileread_helper(filename, 915, [])


def get_experimental_codes():
    print("Exp Codes")


def mio_fileread_helper(filename, root_code, cat_list):
    # print("MIO Helper")

    # code_list = set()
    buffer = ""
    in_term = False
    new_terms = False
    new_term = True

    if not cat_list:
        cat_list = [root_code]

    mio = open(filename, "r")
    for line in mio:
        if in_term and line == "\n":
            in_term = False
            if buffer != "":
                for term_id in cat_list:
                    if buffer.find("is_a: MI:" + ("0" * (4 - len(str(term_id)))) + str(term_id)) != -1:
                        for term_id_2 in cat_list:
                            if ("0" * (4 - len(str(term_id_2))) + str(term_id_2)) == buffer[7:11]:
                                new_term = False
                        if new_term:
                            cat_list.append(int(buffer[7:11]))
                            new_terms = True
                buffer = ""
                new_term = True

        if not in_term and line == "[Term]\n":
            in_term = True

        if in_term:
            if line.startswith("id:") or line.startswith("is_a:"):
                buffer += line

    mio.close()

    print(cat_list)
    print(len(cat_list))

    if new_terms:
        mio_fileread_helper(filename, root_code, cat_list)


def process_biogrid_file():
    print("Process BioGRID")


main()
