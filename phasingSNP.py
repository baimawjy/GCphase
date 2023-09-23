import random
import time


def FM(adjacency_list):
    partition1 = {}
    partition2 = {}
    for vertex in adjacency_list:
        partition1[vertex] = 0
        partition2[vertex] = 1
    iteration = 0
    flag = 1
    snp_list = []
    for vertex in adjacency_list:
        snp_list.append(vertex)

    while flag:
        flag = 0
        n = 0
        iteration += 1
        random.shuffle(snp_list)
        for vertex in snp_list:
            n += 1
            cut_size_old = 0
            cut_size_new = 0
            vertex_direction = partition1[vertex]
            for node in adjacency_list[vertex]:
                node_direction = partition1[node]
                if vertex_direction == node_direction:
                    cut_size_old += adjacency_list[vertex][node][1]
                    cut_size_new += adjacency_list[vertex][node][0]
                else:
                    cut_size_old += adjacency_list[vertex][node][0]
                    cut_size_new += adjacency_list[vertex][node][1]
            if cut_size_new < cut_size_old:
                flag = 1
                partition1[vertex] = partition1[vertex] * -1 + 1
                partition2[vertex] = partition2[vertex] * -1 + 1
    return partition1, partition2


def phasing(output_directory):
    chromosome_list = []
    infile = open(output_directory + "/chromosome_list.txt", 'r')
    for line in infile:
        line = line.strip('\n')
        chromosome_list.append(line)
    infile.close()
    for chromosome in chromosome_list:
        tt = time.localtime(time.time())
        print("chromosome ", chromosome, "start | Now time：", tt[0], "Y", tt[1], "M", tt[2], "D", tt[3], "h", tt[4], "min", tt[5], "s")

        hypergraph = []
        infile = open(output_directory + '/hypergraph/hypergraph-' + chromosome + '.txt', 'r')
        for line in infile:
            line = line.strip('\n')
            edge = line.split('\t')
            hypergraph.append(edge)
        infile.close()
        snp_info = {}
        for edge in hypergraph:
            for node in edge:
                node_position = node[:-2]
                node_direction = node[-1]
                if node_position not in snp_info:
                    if node_direction == '0':
                        snp_info[node_position] = [1, 0]
                    else:
                        snp_info[node_position] = [0, 1]
                else:
                    if node_direction == '0':
                        snp_info[node_position][0] += 1
                    else:
                        snp_info[node_position][1] += 1
        delete_snp = []
        ratio = 0.85
        for snp in snp_info:
            if max(snp_info[snp][0], snp_info[snp][1]) / (snp_info[snp][0] + snp_info[snp][1]) >= ratio or min(
                    snp_info[snp][0], snp_info[snp][1]) == 1:
                delete_snp.append(snp)
        for snp in delete_snp:
            del snp_info[snp]

        hypergraph_clear = []
        for edge in hypergraph:
            edge_clear = []
            for node in edge:
                if node[:-2] in snp_info:
                    edge_clear.append(node)
            hypergraph_clear.append(edge_clear)
        hypergraph = hypergraph_clear
        adjacency_list = {}
        for edge in hypergraph:
            for i in range(len(edge) - 1):
                vertex_position = edge[i][:-2]
                vertex_direction = edge[i][-1]
                if vertex_position not in adjacency_list:
                    adjacency_list[vertex_position] = {}
                for j in range(i + 1, len(edge)):
                    node_position = edge[j][:-2]
                    node_direction = edge[j][-1]
                    if node_position not in adjacency_list[vertex_position]:
                        if node_direction == vertex_direction:
                            adjacency_list[vertex_position][node_position] = [1, 0]
                        else:
                            adjacency_list[vertex_position][node_position] = [0, 1]
                    else:
                        if node_direction == vertex_direction:
                            adjacency_list[vertex_position][node_position][0] += 1
                        else:
                            adjacency_list[vertex_position][node_position][1] += 1
                    if node_position not in adjacency_list:
                        adjacency_list[node_position] = {}
                        if node_direction == vertex_direction:
                            adjacency_list[node_position][vertex_position] = [1, 0]
                        else:
                            adjacency_list[node_position][vertex_position] = [0, 1]
                    else:
                        if vertex_position not in adjacency_list[node_position]:
                            if node_direction == vertex_direction:
                                adjacency_list[node_position][vertex_position] = [1, 0]
                            else:
                                adjacency_list[node_position][vertex_position] = [0, 1]
                        else:
                            if node_direction == vertex_direction:
                                adjacency_list[node_position][vertex_position][0] += 1
                            else:
                                adjacency_list[node_position][vertex_position][1] += 1

        min_edge_len = 1
        delete_vertex = []
        for vertex in adjacency_list:
            delete_node = []
            for node in adjacency_list[vertex]:
                if max(adjacency_list[vertex][node]) - min(adjacency_list[vertex][node]) <= min_edge_len:
                    delete_node.append(node)
            for node in delete_node:
                del adjacency_list[vertex][node]
            if len(adjacency_list[vertex]) == 0:
                delete_vertex.append(vertex)
        for vertex in delete_vertex:
            del adjacency_list[vertex]

        partition1, partition2 = FM(adjacency_list)

        partition = []
        for node in partition1:
            # partition.append(node + ',' + str(partition1[node]))
            partition.append(node + ',0')
        partition = sorted(partition, key=lambda x: int(x.split(',')[0]))

        block_all = [[partition[0]]]
        dif = 2
        error_snp = []
        for i in range(1, len(partition)):
            front_position = block_all[-1][-1][:-2]
            front_direction = str(partition1[front_position])
            next_position = partition[i][:-2]
            next_direction = str(partition1[next_position])
            if next_position not in adjacency_list[front_position]:
                block_all[-1].append(next_position + ',' + next_direction)
                continue
            if front_direction == next_direction:
                if adjacency_list[front_position][next_position][1] - adjacency_list[front_position][next_position][0] >= dif:
                    block_all.append([next_position + ',' + next_direction])
                    error_snp.append([front_position + ',' + front_direction, next_position + ',' + next_direction])
                else:
                    block_all[-1].append(next_position + ',' + next_direction)
            else:
                if adjacency_list[front_position][next_position][0] - adjacency_list[front_position][next_position][1] >= dif:
                    block_all.append([next_position + ',' + next_direction])
                    error_snp.append([front_position + ',' + front_direction, next_position + ',' + next_direction])
                else:
                    block_all[-1].append(next_position + ',' + next_direction)

        for i in range(len(block_all)):
            if i % 2 == 1:
                for node in block_all[i]:
                    node_position = node[:-2]
                    partition1[node_position] = partition1[node_position] * -1 + 1

        block_all = []
        node_flag = {}
        while len(node_flag) != len(adjacency_list):
            queue = []
            block = []
            for vertex in adjacency_list:
                if vertex not in node_flag:
                    node_flag[vertex] = 1
                    block.append(vertex + ',' + str(partition1[vertex]))
                    queue.append(vertex)
                    break
            while len(queue) != 0:
                vertex = queue.pop(0)
                for node in adjacency_list[vertex]:
                    if node not in node_flag:
                        block.append(node + ',' + str(partition1[node]))
                        queue.append(node)
                        node_flag[node] = 1
            block = sorted(block, key=lambda x: int(x[:-2]))
            block_all.append(block)
        block_all = sorted(block_all, key=lambda x: int(x[0][:-2]))

        partition = {}
        snp_list = []
        for block in block_all:
            for node in block:
                snp_list.append(node[:-2])
                partition[node[:-2]] = node[-1]
        snp_list = sorted(snp_list, key=lambda x: int(x))

        block_temp = [[snp_list[0]]]
        dif = 2
        i = 0
        while i < len(snp_list) - 1:
            front_position = snp_list[i]
            front_direction = partition[front_position]
            if i + 1 >= len(snp_list):
                break
            next_position = snp_list[i + 1]
            next_direction = partition[next_position]
            if next_position in adjacency_list[front_position]:
                if front_direction == next_direction:
                    if adjacency_list[front_position][next_position][1] - adjacency_list[front_position][next_position][0] >= dif:
                        block_temp.append([next_position])
                    else:
                        block_temp[-1].append(next_position)
                else:
                    if adjacency_list[front_position][next_position][0] - adjacency_list[front_position][next_position][1] >= dif:
                        block_temp.append([next_position])
                    else:
                        block_temp[-1].append(next_position)
                i += 1
            else:
                if i + 2 >= len(snp_list):
                    block_temp[-1].append(next_position)
                    break
                next_next_position = snp_list[i + 2]
                next_next_direction = partition[next_next_position]
                if next_next_position not in adjacency_list[front_position]:
                    block_temp[-1].append(next_position)
                    i += 1
                    continue

                front_same = 0
                next_same = 0
                for snp_position in adjacency_list[next_position]:
                    snp_direction = partition[snp_position]
                    if int(snp_position) < int(next_position):
                        if next_direction == snp_direction:
                            front_same += adjacency_list[next_position][snp_position][0]
                        else:
                            front_same += adjacency_list[next_position][snp_position][1]
                    else:
                        if next_direction == snp_direction:
                            next_same += adjacency_list[next_position][snp_position][0]
                        else:
                            next_same += adjacency_list[next_position][snp_position][1]

                if front_direction == next_next_direction:
                    if adjacency_list[front_position][next_next_position][1] - \
                            adjacency_list[front_position][next_next_position][0] >= dif:
                        if front_same >= next_same:
                            block_temp[-1].append(next_position)
                            block_temp.append([next_next_position])
                        else:
                            block_temp.append([next_position, next_next_position])
                    else:
                        block_temp[-1].append(next_position)
                        block_temp[-1].append(next_next_position)
                else:
                    if adjacency_list[front_position][next_next_position][0] - \
                            adjacency_list[front_position][next_next_position][1] >= dif:
                        if front_same >= next_same:
                            block_temp[-1].append(next_position)
                            block_temp.append([next_next_position])
                        else:
                            block_temp.append([next_position, next_next_position])
                    else:
                        block_temp[-1].append(next_position)
                        block_temp[-1].append(next_next_position)
                i += 2
        for i in range(len(block_temp)):
            if i % 2 == 1:
                for position in block_temp[i]:
                    partition[position] = str(int(partition[position]) * -1 + 1)

        dif = 2
        block_temp = {}
        for i in range(1, len(snp_list)):
            position = snp_list[i]
            direction = partition[position]
            front_same = 0
            front_different = 0
            front_0 = 0
            front_1 = 0
            for node in adjacency_list[position]:
                if int(node) < int(position):
                    if partition[node] == direction:
                        front_0 += adjacency_list[position][node][0]
                        front_1 += adjacency_list[position][node][1]
                        if adjacency_list[position][node][0] - adjacency_list[position][node][1] >= dif:
                            front_same += 1
                        else:
                            front_different += 1
                    elif partition[node] != direction:
                        front_0 += adjacency_list[position][node][1]
                        front_1 += adjacency_list[position][node][0]
                        if adjacency_list[position][node][1] - adjacency_list[position][node][0] >= dif:
                            front_same += 1
                        else:
                            front_different += 1
            if front_same == 0 and front_different != 0 and front_0 == 0:
                block_temp[position] = 0

        block_temp_temp = [[]]
        for snp in snp_list:
            if snp not in block_temp:
                block_temp_temp[-1].append(snp)
            else:
                block_temp_temp.append([snp])

        for i in range(len(block_temp_temp)):
            if i % 2 == 1:
                for position in block_temp_temp[i]:
                    partition[position] = str(int(partition[position]) * -1 + 1)

        block_all_temp = []
        for block in block_all:
            block_temp = []
            for node in block:
                block_temp.append(node[:-2] + ',' + partition[node[:-2]])
            block_all_temp.append(block_temp)

        outfile = open(output_directory + "/block-FM/block-" + chromosome + ".txt", 'w+')
        for block in block_all_temp:
            temp = '\t'.join(block)
            outfile.write(temp + '\n')
        outfile.close()
        tt = time.localtime(time.time())
        print("chromosome ", chromosome, "done | Now time：", tt[0], "Y", tt[1], "M", tt[2], "D", tt[3], "h", tt[4], "min", tt[5], "s")


def run(output_directory):
    phasing(output_directory)
