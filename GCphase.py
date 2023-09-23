import sys
import os
import get_info
import phasingSNP
import outputResult


if __name__ == '__main__':
    if len(sys.argv) < 7:
        print("Wrong number of parameters")
        sys.exit(1)
    arg = sys.argv
    vcf_file = arg[arg.index('-vcf') + 1]
    bam_file = arg[arg.index('-bam') + 1]
    output_directory = arg[arg.index('-outpath') + 1]

    if not os.path.exists(output_directory + "/hypergraph"):
        os.makedirs(output_directory + "/hypergraph")
    if not os.path.exists(output_directory + "/snpPosition"):
        os.makedirs(output_directory + "/snpPosition")
    if not os.path.exists(output_directory + "/block-FM"):
        os.makedirs(output_directory + "/block-FM")
    if not os.path.exists(output_directory + "/result"):
        os.makedirs(output_directory + "/result")

    get_info.run(vcf_file, bam_file, output_directory)
    phasingSNP.run(output_directory)
    outputResult.run(vcf_file, output_directory)
