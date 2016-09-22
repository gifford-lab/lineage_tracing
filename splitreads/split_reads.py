from collections import namedtuple
from collections import Counter

import argparse
import itertools
import json
import os
import sys
import subprocess
from os  import makedirs
from os.path import exists,join

import Levenshtein as lv

FQRecord = namedtuple("FQRecord", ['name1', 'seq', 'name2', 'qual'])

def get_git_hash():
    script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
    git_hash = subprocess.check_output("cd {0}; git rev-parse --short HEAD".format(script_dir), shell = True).strip()
    return git_hash

def check_index_distance(indexes):
    print >> sys.stderr, "# Check hamming distance for indexes:"
    for index1 in indexes:
        print >> sys.stderr, index1,
        for index2 in indexes:
            print >> sys.stderr, lv.hamming(index1, index2),
        print >> sys.stderr


def fq_iterator(path):
    with open(path) as fq:
        while True:
            name1 = fq.readline().strip()
            if not name1:
                break
            seq = fq.readline().strip()
            name2 = fq.readline().strip()
            qual = fq.readline().strip()
            record = FQRecord(name1, seq, name2, qual)
            yield record

def get_index(record):
    return record.name1.split(':')[-1]

def get_closest_index(query, indexes):
    distances = [(index, lv.hamming(query, index)) for index in indexes]
    distances.sort(key = lambda x: x[1])
    return distances[0] if distances[0][1] != distances[1][1] else (-1,distances[0][1])

def make_out_handles(out_dir, indexes):
    out_handles = {}
    for index in indexes:
        out_handles[index] = {}
        for fq in ['fq1', 'fq2']:
            path = os.path.join(out_dir, '{0}_{1}.fq'.format(index, fq))
            handle = open(path, 'w')
            out_handles[index][fq]= handle
    return out_handles

def close_out_handles(out_handles, indexes):
    for index in indexes:
        for fq in ['fq1', 'fq2']:
            out_handles[index][fq].close()

def write_fq_record(out_handle, fq_record):
    out_handle.write('{0}\n{1}\n{2}\n{3}\n'.format(
        fq_record.name1,
        fq_record.seq,
        fq_record.name2,
        fq_record.qual))

def write_fq_records(out_handles, index, fq1_record, fq2_record):
    write_fq_record(out_handles[index]['fq1'], fq1_record)
    write_fq_record(out_handles[index]['fq2'], fq2_record)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fq1')
    parser.add_argument('-g', '--fq2')
    parser.add_argument('-i', '--indexes')
    parser.add_argument('-o', '--out')
    parser.add_argument('-c', '--hardcutoff',type=int)
    parser.add_argument('-d', '--dryrun', action = "store_true")
    args = parser.parse_args()

    if not exists(args.out):
        makedirs(args.out)

    with open(args.indexes) as f:
        indexes = [line.strip().split()[0] for line in f]
    check_index_distance(indexes)

    fq1_it = fq_iterator(args.fq1)
    fq2_it = fq_iterator(args.fq2)

    counter = Counter()
    read_counter = {index:Counter() for index in indexes}
    i = 0

    out_handles = make_out_handles(args.out, indexes)
    unmap_handle = {'fq1':open(join(args.out,'unmapped_1.fq'),'w'),'fq2':open(join(args.out,'unmapped_2.fq'),'w')}

    for fq1_record, fq2_record in itertools.izip(fq1_it, fq2_it):
        i += 1
        if i % 1000000 == 0:
            print >> sys.stderr, "...Processed {0} reads".format(i)
            if args.dryrun:
                break
        #print fq1_record
        #print fq2_record
        index1 = get_index(fq1_record)
        index2 = get_index(fq2_record)
        assert index1 == index2

        counter[index1] += 1

        closest_index, distance = get_closest_index(index1, indexes)
        if closest_index == -1 or distance > args.hardcutoff:
            write_fq_record(unmap_handle['fq1'],fq1_record)
            write_fq_record(unmap_handle['fq2'],fq2_record)
        else:
            read_counter[closest_index][distance] += 1
            write_fq_records(out_handles, closest_index, fq1_record, fq2_record)

    print >> sys.stderr, '# Number of reads for each index (Hamming distance <= 2):'
    for index in indexes:
        print >> sys.stderr, index,
        s = 0
        for distance in 0, 1, 2:
            print >> sys.stderr, '{0}:{1}'.format(distance, read_counter[index][distance]),
            s += read_counter[index][distance]
        print >> sys.stderr, 'total:{0}'.format(s)

    n_kept_reads = sum([sum([rci for rci in rc.values()]) for
            rc in read_counter.values()])
    print '# Number of reads kept:', n_kept_reads, n_kept_reads/float(i), '%'

    print >> sys.stderr, '# Top indexes observed (exact):'
    for index, count in counter.most_common(10):
        print >> sys.stderr, index, count

    close_out_handles(out_handles, indexes)

    log_path = os.path.join(args.out, 'split_reads.log')
    log = {}
    log['paths'] = {}
    log['paths']['fq1'] = args.fq1
    log['paths']['fq2'] = args.fq2
    log['paths']['indexes'] = args.indexes
    log['read_counter'] = read_counter
    log['git_hash'] = get_git_hash()
    with open(log_path, 'w') as f:
        json.dump(log, f)

if __name__ == '__main__':
    main()
