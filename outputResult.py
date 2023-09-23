import time


def modify(vcf_file, output_directory):
    tt = time.localtime(time.time())
    print("output file start | Now time：", tt[0], "Y", tt[1], "M", tt[2], "D", tt[3], "h", tt[4], "min", tt[5], "s")
    snp_all = {}
    chromosome_list = []
    infile = open(output_directory + "/chromosome_list.txt", 'r')
    for line in infile:
        line = line.strip('\n')
        chromosome_list.append(line)
        snp_all[line] = {}
    infile.close()
    for chromosome in chromosome_list:
        infile = open(output_directory + '/block-FM/block-' + chromosome + '.txt', 'r')
        for line in infile:
            line = line.strip('\n').split('\t')
            block_first_snp = line[0][:-2]
            for snp in line:
                snp_position = snp[:-2]
                snp_direction = snp[-1]
                snp_all[chromosome][snp_position] = [snp_direction, block_first_snp]
        infile.close()

    infile = open(vcf_file, 'r')
    outfile = open(output_directory + '/result/mine.vcf', 'w+')
    for line in infile:
        line = line.strip('\n')
        if line[0] == '#':
            outfile.write(line + '\n')
        else:
            line = line.split('\t')
            line_snp_position = line[1]
            line_snp_chromosome = line[0]
            if line_snp_position in snp_all[line_snp_chromosome]:
                line_snp_info = line[-1].split(':')
                line_snp_flag = line[-2].split(':')
                snp_phase_num = 0
                snp_block_num = -1
                for i in range(len(line_snp_flag)):
                    if line_snp_flag[i] == 'GT':
                        snp_phase_num = i
                    elif line_snp_flag[i] == 'PS':
                        snp_block_num = i
                if snp_all[line_snp_chromosome][line_snp_position][0] == '0':
                    line_snp_info[snp_phase_num] = '0|1'
                elif snp_all[line_snp_chromosome][line_snp_position][0] == '1':
                    line_snp_info[snp_phase_num] = '1|0'
                line_snp_info[snp_block_num] = snp_all[line_snp_chromosome][line_snp_position][1]
                temp = ':'.join(line_snp_info)
                line[-1] = temp
                temp = '\t'.join(line)
                outfile.write(temp + '\n')
            else:
                temp = '\t'.join(line)
                outfile.write(temp + '\n')
    infile.close()
    outfile.close()
    tt = time.localtime(time.time())
    print("output file done | Now time：", tt[0], "Y", tt[1], "M", tt[2], "D", tt[3], "h", tt[4], "min", tt[5], "s")


def run(vcf_file, output_directory):
    modify(vcf_file, output_directory)