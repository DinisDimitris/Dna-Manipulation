# this file is not part of the assessment and serves only for the purpose of creating data to time the scripts

from random import randint

def generate_dna_sequence(l):
    alphabet = "CTAG"

    dna_sequence = ""

    for i in range(l):
        dna_sequence += alphabet[randint(0, 3)]

    return dna_sequence


def find_restriction_enzyme(dna_seq, l, find_after, num_enz):

    position = randint(find_after, find_after-enzyme_length+len(dna_seq)//num_enz)

    restriction_enzyme = dna_seq[position:position+l]

    return restriction_enzyme, position + enzyme_length


with open("timing_test_data.txt", "w") as f_out:
    for i in range(1):

        dna_length = 10000000
        enzyme_length = 10
        strand_length = 10
        number_of_enzymes = 5

        dna = generate_dna_sequence(dna_length)
        f_out.write("%d %d\n%s\n" % (dna_length, number_of_enzymes, dna))
        last_enzyme_match_position = 0

        for k in range(number_of_enzymes):
            enzyme, last_enzyme_match_position = find_restriction_enzyme(dna, enzyme_length, last_enzyme_match_position, number_of_enzymes)
            splicing_index = randint(0, enzyme_length)
            strand = generate_dna_sequence(strand_length)

            f_out.write("%d %d %s %d %s\n" % (enzyme_length, strand_length, enzyme, splicing_index, strand))

        f_out.write("\n\n")
