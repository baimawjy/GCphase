import pysam
import time


def get_snp_position(vcf_file_path, output_directory):
    chromosome_list = []
    vcf_file = open(vcf_file_path, 'r')
    chromosome = ""
    SNP_position = []
    SNP_chr_position = []
    line_split = []
    # tt = time.time()
    tt = time.localtime(time.time())
    print("get SNP info start | Now time:", tt[0], "Y", tt[1], "M", tt[2], "D", tt[3], "h", tt[4], "min", tt[5], "s")
    for line in vcf_file:
        if line[0] != '#':
            line_split = line.split('\t')
            if line_split[-1][0] == '1':
                continue
            if line_split[0] != chromosome:
                if chromosome == "":
                    chromosome = line_split[0]
                else:
                    SNP_position.append(SNP_chr_position)
                    chromosome = line_split[0]
                    SNP_chr_position = []
            SNP_chr_position.append([line_split[0], line_split[1], line_split[3], line_split[4]])
    if len(SNP_chr_position) != 0:
        SNP_position.append(SNP_chr_position)
        SNP_chr_position.append([line_split[0], line_split[1], line_split[3], line_split[4]])
    tt = time.localtime(time.time())
    print("get SNP info done | Now time", tt[0], "Y", tt[1], "M", tt[2], "D", tt[3], "h", tt[4], "min", tt[5], "s")
    print("chromosome num:", len(SNP_position))
    for snp_chromosome in SNP_position:
        chromosome = snp_chromosome[0][0]
        chromosome_list.append(chromosome)
        outfile = open(output_directory + "/snpPosition/snp-" + chromosome + ".txt", 'w+')
        for snp in snp_chromosome:
            outfile.write(snp[1] + ',' + snp[2] + ',' + snp[3] + '\n')
        outfile.close()
    outfile = open(output_directory + '/chromosome_list.txt', 'w+')
    for chromosome in chromosome_list:
        outfile.write(chromosome + '\n')
    outfile.close()
    return SNP_position, chromosome_list


def get_read_snp_coverage(SNP_position, bam_file_path, output_directory):
    t = time.time()
    tt = time.localtime(time.time())
    print("get reads info start | Now time", tt[0], "Y", tt[1], "M", tt[2], "D", tt[3], "h", tt[4], "min", tt[5], "s")
    chromosome_flag = {}
    hypergraph_all = []
    for i in range(len(SNP_position)):
        snp_chromosome = SNP_position[i]
        chromosome = snp_chromosome[0][0]
        chromosome_flag[chromosome] = i
        hypergraph_all.append([])

    num_alignments = 0

    for reads in pysam.AlignmentFile(bam_file_path, "rb").fetch(until_eof=True):
        num_alignments += 1
        print('\r', "line num:", num_alignments, end='')
        flag = reads.flag
        if flag != 0 and flag != 16:
            continue
        chromosome = reads.reference_name

        if chromosome not in chromosome_flag:
            continue
        snp_num = len(SNP_position[chromosome_flag[chromosome]])
        read_start = reads.reference_start + 1
        read_end = reads.reference_end + 1
        read_snp_position = []
        read_snp_position_in_ref = []
        left = 0
        right = snp_num - 1
        position = -1
        while left <= right:
            middle = (left + right) // 2
            temp = int(SNP_position[chromosome_flag[chromosome]][middle][1])
            if temp > read_start:
                right = middle - 1
            elif temp < read_start:
                left = middle + 1
            elif temp == read_start:
                position = middle
                break
        if position == -1:
            position = right + 1

        while position < snp_num and int(SNP_position[chromosome_flag[chromosome]][position][1]) <= read_end:
            read_snp_position.append(int(SNP_position[chromosome_flag[chromosome]][position][1]))
            read_snp_position_in_ref.append(position)
            position += 1
        snp_num_read = len(read_snp_position)
        if snp_num_read <= 1:
            continue
        reads_seq = reads.query_sequence
        cigar_tuples = reads.cigartuples
        read_coverage = read_start
        offset = 0
        position_n = 0
        read_snp_base = []
        for op in cigar_tuples:
            cigar_flag = op[0]
            cigar_len = op[1]
            if cigar_flag == 0:  # M
                read_coverage += cigar_len
                while position_n < snp_num_read and read_snp_position[position_n] < read_coverage:
                    read_snp_base.append(reads_seq[read_snp_position[position_n] - read_start + offset])
                    position_n += 1
                if position_n == snp_num_read:
                    break
            elif cigar_flag == 1:  # I
                offset += cigar_len
            elif cigar_flag == 2:  # D
                read_coverage += cigar_len
                while position_n < snp_num_read and read_snp_position[position_n] < read_coverage:
                    read_snp_base.append('-')
                    position_n += 1
                if position_n == snp_num_read:
                    break
                offset -= cigar_len
            elif cigar_flag == 4:  # S
                offset += cigar_len
        edge = []
        for i in range(len(read_snp_base)):
            read_base = read_snp_base[i]
            snp_info = SNP_position[chromosome_flag[chromosome]][read_snp_position_in_ref[i]]
            if read_base == snp_info[2]:
                edge.append(snp_info[1] + ',0')
            elif read_base == snp_info[3]:
                edge.append(snp_info[1] + ',1')
        if len(edge) <= 1:
            continue
        hypergraph_all[chromosome_flag[chromosome]].append(edge)

    for chromosome in chromosome_flag:
        outfile = open(output_directory + "/hypergraph/hypergraph-" + chromosome + ".txt", 'w+')
        hypergraph = hypergraph_all[chromosome_flag[chromosome]]
        for edge in hypergraph:
            temp = '\t'.join(edge)
            outfile.write(temp + '\n')
        outfile.close()

    tt = time.localtime(time.time())
    print("get reads info done | Now time：", tt[0], "Y", tt[1], "M", tt[2], "D", tt[3], "h", tt[4], "min", tt[5], "s")


def run(vcf_file, bam_file, output_directory):
    SNP_position, chromosome_list = get_snp_position(vcf_file, output_directory)
    get_read_snp_coverage(SNP_position, bam_file, output_directory)
