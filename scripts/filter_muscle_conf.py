import copy
from collections import defaultdict
import argparse
import logging

__author__ = "Bernhard Bein, 2025"

## Logging
log = logging.getLogger(__name__)

DESCRIPTION = '''
Script to filter alignments created with muscle, exchanging bases below user-specified confidence score thresholds with gap characters ('-').
Needs a muscle confidence score file and a alignment of the same length as input. 
'''

def argument_parser():
  """Parse CMD args."""
  app = argparse.ArgumentParser(description=DESCRIPTION, formatter_class= lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=6, indent_increment=2))
  
  app.add_argument(
  "-ia", 
  "--input_alignment",
  action="store",
  dest="input",
  type=str,
  help=
'''Input fasta file created by a muscle.
  ''')

  app.add_argument(
  "-ic", 
  "--input_confidence",
  action="store",
  dest="conf",
  type=str,
  help=
'''Input letter confidence file created by muscle.
  ''')

  app.add_argument(
  "-o", 
  "--output_fasta",
  action="store",
  dest="output",
  type=str,
  help=
'''Output name of the cleaned fasta file.
  ''')

  app.add_argument(
  "-t", 
  "--threshold",
  action="store",
  dest="threshold",
  type=str,
  help=
'''confidence threshold below which bases will be turned into gaps '-'. 
For example, if the user provides -t 9, all bases with a letter confidence score below 9 will be removed.
  ''')

  args = app.parse_args()
  return args

def read_fasta(fasta_path):
    fasta_dict = defaultdict(str)
    with open(fasta_path, "r") as fasta:
        for line in fasta:
            fasta_line  = line.strip()
            if fasta_line.startswith(">"):
                header = fasta_line
            else:
                body = fasta_line
                fasta_dict[header] += body
    return fasta_dict

def filter_alignment(fasta, confidence, threshold, reference):
    filtered_fasta = {}
    removed_characters = 0
    for species in fasta.keys():
        if species == ">" + reference:
            log.info("Reference was skipped for base removal")
            filtered_fasta[species] = fasta[species]
            continue
        if confidence[species]:
            current_seq = ''
            for i, nuc in enumerate(fasta[species]):
                if confidence[species][i] == "-" or not isinstance(confidence[species][i], int):
                    current_seq += "-"
                elif int(confidence[species][i]) < threshold:
                    current_seq += "-"
                    removed_characters += 1
                else:
                    current_seq += nuc
            filtered_fasta[species] = current_seq
        log.info("for species {}, {} characters where removed.".format(species, removed_characters))
        removed_characters = 0 
    return filtered_fasta

## CLI args masked for snakemake executrion
def main() -> None:

    input_fasta = snakemake.input.in_ali
    output = snakemake.output.filtered_ali
    input_confidence = snakemake.input.conf_scores
    reference = snakemake.params.reference
    threshold = 9
    logfile = snakemake.log[0]

    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.DEBUG,
        format="[%(levelname)s] %(message)s",
    )

    log.info("The user-set threshold for letter confidence is {}".format(threshold))
    fasta_dict = read_fasta(input_fasta)
    conf_dict = read_fasta(input_confidence)
    filtered_dict = filter_alignment(fasta_dict, conf_dict, threshold, reference)
    with open(output, "w") as out:
        for spec in filtered_dict.keys():
            out.write(spec+"\n"+filtered_dict[spec]+"\n")


if __name__ == "__main__":
    main()
            
