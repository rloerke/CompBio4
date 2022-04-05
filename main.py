"""
Comp Bio Assignment 4
Written by Ray Loerke
BioGRID MITAB File Processor
"""


def main():
    code_list = format_helper("mi.owl")

    phys_set = get_physical_interaction_codes(code_list)
    print(phys_set)
    print(len(phys_set))
    exp_set = get_experimental_codes(code_list)
    print(exp_set)
    print(len(exp_set))


def get_physical_interaction_codes(code_list):
    print("Phys Codes")
    return mio_helper(code_list, 915, set())


def get_experimental_codes(code_list):
    print("Exp Codes")
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


def process_biogrid_file():
    print("Process BioGRID")


main()
