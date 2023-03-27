"""
Calculate minimal distances between sequences in an alignment and a set of focal sequences
"""
import argparse
from augur.io import read_sequences
from random import shuffle
from collections import defaultdict
import Bio
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio import AlignIO, SeqIO
from scipy import sparse
import sys


def compactify_sequences(sparse_matrix, sequence_names):
    sequence_groups = defaultdict(list)
    for s, snps in zip(sequence_names, sparse_matrix):
        ind = snps.nonzero()
        vals = np.array(snps[ind])
        if len(ind[1]):
            sequence_groups[tuple(zip(ind[1], vals[0]))].append(s)
        else:
            sequence_groups[tuple()].append(s)

    return sequence_groups

INITIALISATION_LENGTH = 1000000 # 1Mb的长度

def sequence_to_int_array(s, fill_value=110, fill_gaps=True):
    '''
    把字符串转换为小写字母的ascii值。
    acgt分别为97, 99, 103, 116。-为45，n为110
    fill_gaps为true，则用n填充-
    '''
    seq = np.frombuffer(str(s).lower().encode('utf-8'), dtype=np.int8).copy()
    if fill_gaps:
        seq[(seq!=97) & (seq!=99) & (seq!=103) & (seq!=116)] = fill_value
    else:
        seq[(seq!=97) & (seq!=99) & (seq!=103) & (seq!=116) & (seq!=45)] = fill_value
    return seq

# Function adapted from https://github.com/gtonkinhill/pairsnp-python
def calculate_snp_matrix(fastafile, consensus=None, zipped=False, fill_value=110, chunk_size=0, ignore_seqs=None):
    # This function generate a sparse matrix where differences to the consensus are coded as integers.
    if ignore_seqs is None:
        ignore_seqs = []

    # 以下初始化数组形状，但是值是随机的
    row = np.empty(INITIALISATION_LENGTH)
    col = np.empty(INITIALISATION_LENGTH, dtype=np.int64)
    val = np.empty(INITIALISATION_LENGTH, dtype=np.int8)

    r = 0
    n_snps = 0 # 分析的snps总数目
    nseqs = 0 # 分析的序列数目
    seq_names = [] # 分析的序列名
    filled_positions = []
    current_length = INITIALISATION_LENGTH

    for record in fastafile:
        h = record.name
        s = str(record.seq)

        if h in ignore_seqs:
            continue
        if consensus is None:
            align_length = len(s)
            # Take consensus as first sequence
            consensus = sequence_to_int_array(s, fill_value=fill_value)
        else:
            align_length = len(consensus)

        nseqs +=1
        seq_names.append(h)

        if(len(s)!=align_length):
            raise ValueError('Fasta file appears to have sequences of different lengths!')

        s = sequence_to_int_array(s, fill_value=fill_value)
        snps = (consensus!=s) & (s!=fill_value) # ref与s不一致，并且s的碱基不是n，True or False
        right = n_snps + np.sum(snps)
        filled_positions.append(np.where(s==fill_value)[0]) 
        # np.where的效果等同于np.asarray(s==fill_value).nonzero()[0]，按照维度返回True的位置
        # 在这里，返回s中n所在的索引

        if right >= (current_length/2):
            current_length = current_length + INITIALISATION_LENGTH
            # 以下调整大小，并用0填充
            row.resize(current_length)
            col.resize(current_length)
            val.resize(current_length)

        row[n_snps:right] = r # 第1条序列的位置都为0，第2条序列的位置都为1
        col[n_snps:right] = np.flatnonzero(snps) # 给出所有snps的索引
        val[n_snps:right] = s[snps] # 返回所有snp位点的碱基值
        r += 1
        n_snps = right
        if chunk_size and chunk_size==nseqs:
            # 处理足够多的序列
            break

    if nseqs==0:
        return None

    row = row[0:right]
    col = col[0:right]
    val = val[0:right]

    # 创建稀疏矩阵，值为val数组
    sparse_snps = sparse.csc_matrix((val, (row, col)), shape=(nseqs, align_length))

    return {'snps': sparse_snps, 'consensus': consensus, 'names': seq_names, 'filled_positions': filled_positions}

# Function adapted from https://github.com/gtonkinhill/pairsnp-python
def calculate_distance_matrix(sparse_matrix_A, sparse_matrix_B, consensus):

    n_seqs_A = sparse_matrix_A.shape[0]
    n_seqs_B = sparse_matrix_B.shape[0]

    d = (1*(sparse_matrix_A==97)) * (sparse_matrix_B.transpose()==97) # a
    d = d + (1*(sparse_matrix_A==99) * (sparse_matrix_B.transpose()==99)) # c
    d = d + (1*(sparse_matrix_A==103) * (sparse_matrix_B.transpose()==103)) # g
    d = d + (1*(sparse_matrix_A==116) * (sparse_matrix_B.transpose()==116)) # t

    # (n, len)的矩阵A与(len, m)的矩阵B，相乘->(n, m)的结果矩阵，表示n条序列与m条序列相比，其SNP位点是否相同
    d = d.todense()

    # 计算n碱基是否相同，获得比较矩阵
    n_comp = (1*(sparse_matrix_A==110) * ((sparse_matrix_B==110).transpose())).todense()
    d = d + n_comp
    
    # 计算A或B中SNP和N的数目
    temp_total = np.zeros((n_seqs_A, n_seqs_B))
    temp_total[:] = (1*(sparse_matrix_A>0)).sum(1)
    temp_total += (1*(sparse_matrix_B>0)).sum(1).transpose()

    # 计算A和B共有的SN和N位点数目
    total_differences_shared = (1*(sparse_matrix_A>0)) * (sparse_matrix_B.transpose()>0)

    # 计算A或B中n碱基的数目
    n_total = np.zeros((n_seqs_A, n_seqs_B))
    n_sum = (1*(sparse_matrix_A==110)).sum(1)
    n_total[:] = n_sum
    n_total += (1*(sparse_matrix_B==110)).sum(1).transpose()

    # 计算A或B独有的n
    diff_n = n_total - 2*n_comp
    d = temp_total - total_differences_shared.todense() - d - diff_n

    # 距离计算公式为：A中SNP和N的数目 + B中SNP和N的数目 - AB共有的SNP和N的数目 - AB中SNP和N碱基相同的数目 - A和B中独有的N的数目

    return d

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="generate priorities files based on genetic proximity to focal sample",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", type=str, required=True, help="FASTA file of alignment")
    parser.add_argument("--reference", type = str, required=True, help="reference sequence (FASTA)")
    parser.add_argument("--ignore-seqs", type = str, nargs='+', help="sequences to ignore in distance calculation")
    parser.add_argument("--focal-alignment", type = str, required=True, help="focal sample of sequences")
    parser.add_argument("--chunk-size", type=int, default=10000, help="number of samples in the global alignment to process at once. Reduce this number to reduce memory usage at the cost of increased run-time.")
    parser.add_argument("--output", type=str, required=True, help="FASTA file of output alignment")
    args = parser.parse_args()

    # load entire alignment and the alignment of focal sequences (upper case -- probably not necessary)
    ref = sequence_to_int_array(SeqIO.read(args.reference, 'fasta').seq)
    alignment_length = len(ref)

    # 读取focal fasta文件，返回Bio.SeqRecord.SeqRecord对象
    focal_seqs = read_sequences(args.focal_alignment)
    focal_seqs_dict = calculate_snp_matrix(focal_seqs, consensus = ref, ignore_seqs=args.ignore_seqs)

    if focal_seqs_dict is None:
        print(
            f"ERROR: There are no valid sequences in the focal alignment, '{args.focal_alignment}', to compare against the full alignment.",
            "Check your subsampling settings for the focal alignment or consider disabling proximity-based subsampling.",
            file=sys.stderr
        )
        sys.exit(1)

    seqs = read_sequences(args.alignment)

    # export priorities
    fh_out = open(args.output, 'w')
    fh_out.write('strain\tclosest strain\tdistance\n')

    chunk_size=args.chunk_size
    chunk_count = 0
    while True:
        context_seqs_dict = calculate_snp_matrix(seqs, consensus=ref, chunk_size=chunk_size)
        if context_seqs_dict is None:
            break

        print("Reading the alignments.", chunk_count*chunk_size)

        # calculate number of masked sites in either set
        # 计算序列中N的数目
        mask_count_focal = np.array([len(x) for x in focal_seqs_dict['filled_positions']])
        mask_count_context = {s: len(x) for s,x in zip(context_seqs_dict['names'], context_seqs_dict['filled_positions'])}

        # for each context sequence, calculate minimal distance to focal set, weigh with number of N/- to pick best sequence
        d = np.array(calculate_distance_matrix(context_seqs_dict['snps'], focal_seqs_dict['snps'], consensus = context_seqs_dict['consensus']))
        closest_match = np.argmin(d+mask_count_focal/alignment_length, axis=1) # 对于每个context，选择最近的focal
        print("Done finding closest matches.")

        minimal_distance_to_focal_set = {}
        for context_index, focal_index in enumerate(closest_match):
            # 迭代每条context序列，输出最近的focal序列id及距离
            minimal_distance_to_focal_set[context_seqs_dict['names'][context_index]] = (d[context_index, focal_index], focal_seqs_dict["names"][focal_index])

        for seqid in minimal_distance_to_focal_set:
            fh_out.write(f"{seqid}\t{minimal_distance_to_focal_set[seqid][1]}\t{minimal_distance_to_focal_set[seqid][0]}\n")

        chunk_count += 1

    fh_out.close()
