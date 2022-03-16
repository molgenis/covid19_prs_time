#!/usr/env/bin python3

import sys
import os
import argparse
import gzip

__author__ = "C.A. (Robert) Warmerdam"
__credits__ = ["Robert Warmerdam"]
__version__ = "0.0.2"
__status__ = "Development"


class ArgumentParser:

    def __init__(self):
        self.parser = self.create_argument_parser()
        self.add_column_indices_arguments()
        self.add_output_file_path_prefix_argument()
        self.add_snp_db_file_path_argument()
        self.add_summary_statistics_file_path_argument()

    def parse_input(self, argv):
        """
        Parse command line input.
        :param argv: given arguments
        :return: parsed arguments
        """
        args = self.parser.parse_args(argv)
        return args

    def create_argument_parser(self):
        """
        Method creating an argument parser
        :return: parser
        """
        parser = argparse.ArgumentParser(description="Summary stats annotation",
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
        return parser

    def add_summary_statistics_file_path_argument(self):
        self.parser.add_argument('-s', '--ssf-path', type=self.is_readable_file,
                                 required=True,
                                 default=None,
                                 help="Path pointing to summary statistics file")

    def add_snp_db_file_path_argument(self):
        self.parser.add_argument('-d', '--snp-db-file-path', type=self.is_readable_file,
                                 required=True,
                                 default=None,
                                 help="File of variants with 5 columns representing the "
                                      "chromosome, bp, first allele, alternative allele and identifier.")

    def add_output_file_path_prefix_argument(self):
        self.parser.add_argument('-o', '--output-file', type=argparse.FileType('wb'),
                                 required=False, default=None,
                                 help="Output directory")

    def add_column_indices_arguments(self):
        self.parser.add_argument('-CHR', '--chromosome', type=int, required=True, default=None)
        self.parser.add_argument('-BP', '--position', type=int, required=True, default=None)
        self.parser.add_argument('-A1', '--effect-allele', type=int, required=True, default=None)
        self.parser.add_argument('-A2', '--non-effect-allele', type=int, required=True, default=None)

    @staticmethod
    def is_boolean_value(value):
        if isinstance(value, bool):
            return value
        if value.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif value.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    @staticmethod
    def is_readable_file(file_path):
        """
        Checks whether the given directory is readable
        :param file_path: a path to a file in string format
        :return: file_path
        :raises: Exception: if the given path is invalid
        :raises: Exception: if the given directory is not accessible
        """
        if not os.path.isfile(file_path):
            raise Exception("file path:{0} is not a valid file path".format(file_path))
        if os.access(file_path, os.R_OK):
            return file_path
        else:
            raise Exception("file path:{0} is not a readable file".format(file_path))


class FlippedAnnotator():
    def __init__(self, snp_db_file_path):
        self.snp_db_file_path = snp_db_file_path

    def find_annotations(self, queries, annotations):
        #print(annotations, file=sys.stderr)
        #print(queries, file=sys.stderr)
        with gzip.open(self.snp_db_file_path, 'rt') as opened_snp_db:
            number_of_variants = 0
            processed_queries = 0
            chromosomes = dict()
            print("Opened '{}'...".format(self.snp_db_file_path), file=sys.stderr)

            for line in opened_snp_db:
                # Only process the line if the lien starts with '##'
                if not line.startswith("#"):
                    #print(line, file=sys.stderr)
                    split_line = line.split()
                    #print(split_line, file=sys.stderr)

                    if not split_line[0] in chromosomes:
                        chrom = split_line[0].replace('_', ' ').replace('.', ' ').split()[1].lstrip("0")
                        chromosomes[split_line[0]] = chrom

                    chrom = chromosomes[split_line[0]]
                    basepair = split_line[1]
                    allele_1 = split_line[3]
                    alleles_alt = split_line[4]

                    # snp_db[":".join(split_line[:4])] = split_line[4]
                    for allele_2 in alleles_alt.split(","):
                        variant = ":".join([chrom,
                                         basepair,
                                         allele_1,
                                         allele_2])

                        #print(variant, file=sys.stderr)

                        if variant in queries:
                            annotations[queries[variant]] = split_line[2]
                            processed_queries += 1

                    number_of_variants += 1

                    if number_of_variants % 1000000 == 0:
                        print("Number of processed variants: {}, processed queries {}"
                              .format(number_of_variants, processed_queries), file=sys.stderr)
                        # for key, value in snp_db.items():
                        #    print(key, value, file=sys.stderr)
        return annotations


# class Annotator:
#
#     def __init__(self, snp_db):
#         self.snp_db = snp_db
#
#     @classmethod
#     def from_snp_db_file(cls, snp_db_file_path):
#         snp_db = dict()
#
#         with open(snp_db_file_path) as opened_snp_db:
#             number_of_variants = 0
#             chromosomes = dict()
#             print("Opened '{}'...".format(snp_db_file_path), file=sys.stderr)
#
#             for line in opened_snp_db:
#                 # Only process the line if the lien starts with '##'
#                 if not line.startswith("#"):
#                     split_line = line.split()
#
#                     if not split_line[0] in chromosomes:
#                         chrom = split_line[0].replace('_', ' ').replace('.', ' ').split()[1].lstrip("0")
#                         chromosomes[split_line[0]] = chrom
#
#                     chrom = chromosomes[split_line[0]]
#                     basepair = split_line[1]
#                     allele_1 = split_line[3]
#                     alleles_alt = split_line[4]
#
#                     # snp_db[":".join(split_line[:4])] = split_line[4]
#                     for allele_2 in alleles_alt.split(","):
#                         snp_db[":".join([chrom,
#                                          basepair,
#                                          allele_1,
#                                          allele_2])] = split_line[2]
#                         number_of_variants += 1
#
#                     if number_of_variants % 1000000 == 0:
#                         print("Number of processed variants:", number_of_variants, file=sys.stderr)
#                         # for key, value in snp_db.items():
#                         #    print(key, value, file=sys.stderr)
#         return cls(snp_db)


def get_value(line, index):
    if index < 0:
        return "-"
    else:
        return line[index]


def main(argv=None):
    if argv is None:
        argv = sys.argv

    #print("WHOOP")

    args = ArgumentParser().parse_input(argv[1:])

    print("Starting with arguments...", file=sys.stderr)
    print(args, file=sys.stderr)

    # snp_db_annotator = Annotator.from_snp_db_file(args.snp_db_file_path)
    snp_db_annotator = FlippedAnnotator(args.snp_db_file_path)

    print("Created annotator...", file=sys.stderr)

    not_found = 0

    effect_allele_index = args.effect_allele - 1
    non_effect_allele_index = args.non_effect_allele - 1
    chromosome_index = args.chromosome - 1
    position_index = args.position - 1

    queries = dict()
    annotations = list()

    with open(args.ssf_path) as opened_summary:
        #print(opened_summary.readline().rstrip() + "\tRSID")
        index = 0

        for line in opened_summary:
            split_line = line.strip().split("\t")

            effect_allele = split_line[effect_allele_index].upper()
            non_effect_allele = split_line[non_effect_allele_index].upper()
            
            key_frst = ":".join([split_line[chromosome_index],
                                 split_line[position_index],
                                 effect_allele,
                                 non_effect_allele])

            key_scnd = ":".join([split_line[chromosome_index],
                                 split_line[position_index],
                                 non_effect_allele,
                                 effect_allele])

            # identifier = snp_db_annotator.snp_db.get(
            #     key_frst, snp_db_annotator.snp_db.get(
            #         key_scnd, split_line[4]))

            queries[key_frst] = index
            queries[key_scnd] = index

            index += 1
            annotations.append("_".join([":".join([split_line[chromosome_index],
                                 split_line[position_index]]),
                                 non_effect_allele,
                                 effect_allele]))

            # if split_line[4] == "12:53118797:T:C":
            #    print("12:53118797:T:C", key_frst, key_scnd, identifier, file=sys.stderr)

            # print(split_snp_label)
            # print("\t".join([line.rstrip(), identifier]))

    annotations = snp_db_annotator.find_annotations(queries, annotations)
    
    #print("RSID")
    
    for annotation in annotations:
        print(annotation)
    return 0


if __name__ == "__main__":
    sys.exit(main())
