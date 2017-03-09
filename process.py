from os import listdir
import os
import sys
from subprocess import call
from os.path import expanduser

UNZIPPED_DIRECTORY = "1_unzipped/"
JOINED_DIRECTORY = "2_joined/"
CONCAT_DIRECTORY = "3_concat_filter/"
TRIMMED_DIRECTORY = "4_trimmed/"
UNIQUES_DIRECTORY = "5_uniques/"
OTU_DIRECTORY = "6_1_OTUs/"
OTU_PARSE_DIRECTORY = "6_2_uparseout/"
TABLE_DIRECTORY = "7_tables/"

folders = [UNZIPPED_DIRECTORY, JOINED_DIRECTORY, CONCAT_DIRECTORY, UNIQUES_DIRECTORY, TRIMMED_DIRECTORY, UNIQUES_DIRECTORY, OTU_DIRECTORY, OTU_PARSE_DIRECTORY, TABLE_DIRECTORY]

F_PRIMER = 'GTGCCAGCMGCCGCGGTAA'
F_OFFSET = 3
R_PRIMER = 'GGACTACHVGGGTWTCTAAT'
MIN_LENGTH_THRESHOLD = 90

walk_dir = sys.argv[1]

def makeFolders(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
makeFolders(folders)

prefix = 'ssr_'


def unzip():
    for root, subdirs, files in os.walk(walk_dir):
        if len(files) > 0:
            for f in files:
                if f.split('.')[-1] == 'gz':
                    print 'Copying/Unzipping ' + f
                    call(["cp", os.path.join(root, f), UNZIPPED_DIRECTORY])
                    identifier = f[9:]
                    newname = os.path.join(root, f).split('/')[-2] + identifier
                    call(["mv", UNZIPPED_DIRECTORY + f, UNZIPPED_DIRECTORY + newname])
                    call(["gunzip", UNZIPPED_DIRECTORY + newname, '-f'])

def pre_join():
    fastqFiles = []
    for f in listdir(UNZIPPED_DIRECTORY):
        if f.split('.')[len(f.split('.'))-1] == 'fastq':
            fastqFiles.append(f)

    forwardReads = []
    for f in fastqFiles:
        if f.split('__')[1] == 'R1':
            forwardReads.append(f)
    forwardReads = sorted(forwardReads)

    pairedReads = []
    for fwd in forwardReads:
        ID = fwd.split('__')[:1] + [fwd.split('__')[-1]]
        for rev in fastqFiles:
            if rev.split('__')[:1]  + [rev.split('__')[-1]] == ID and rev.split('__')[1] == 'R2':
                pairedReads.append([fwd, rev])
    print 'Found ' + str(len(pairedReads)) + ' pairs'
    for i in range(len(pairedReads)):
        print 'joining ' + str(i) + '/' + str(len(pairedReads)) + ' :' + str(pairedReads[i])
        join(pairedReads[i])

def join(pair):
    call(["usearch", '-fastq_join',
        UNZIPPED_DIRECTORY + pair[0],
        '-reverse', UNZIPPED_DIRECTORY + pair[1],
        '-fastqout',
        JOINED_DIRECTORY + pair[0].split('__R1__')[0] + '_R1R2_' + pair[0].split('__R1__')[-1].split('.')[0] + 'joined.fastq'])


def L00_concat_filter():
    prefixes = set()
    joined = listdir(JOINED_DIRECTORY)

    trueJoined = []
    for f in joined:
        if 'joined' in f:
            trueJoined.append(f)
    joined = trueJoined

    for f in joined:
        prefixes.add(f.split('L00')[0])

    for pre in prefixes:
        cats = []
        for f in joined:
            if pre in f:
                cats.append(JOINED_DIRECTORY + f)

        command = 'cat '
        for c in cats:
            command += c
            command += ' '
        command += '> ' + CONCAT_DIRECTORY + pre + 'L001.concat.fastq'
        print 'concatenating', cats
        call(command, shell=True)


        command = 'usearch -fastq_filter ' + CONCAT_DIRECTORY + pre + 'L001.concat.fastq -relabel cluster. -fastq_minlen ' + str(MIN_LENGTH_THRESHOLD) + ' -fastqout ' + CONCAT_DIRECTORY + pre + 'L001.filtered.fastq'
        print 'filtering', cats
        call(command, shell=True)


def trim():
    filtered = []
    for f in listdir(CONCAT_DIRECTORY):
        if 'filtered' in f:
            filtered.append(f)

    for f in filtered:
        print filtered.index(f), len(filtered)
        call(['usearch',
            '-fastx_truncate', CONCAT_DIRECTORY + f,
            '-stripleft', str(len(F_PRIMER) + F_OFFSET),
            '-stripright', str(len(R_PRIMER) + 4),
            '-fastaout', TRIMMED_DIRECTORY + f.split('.')[0] + '.trimmed.fasta'
            ])

def unique():
    trimmed = []
    for f in listdir(TRIMMED_DIRECTORY):
        if 'trimmed' in f:
            trimmed.append(f)

    for f in trimmed:
        print trimmed.index(f), len(trimmed)
        call(['usearch',
            '-fastx_uniques', TRIMMED_DIRECTORY + f,
            '-sizeout', '-relabel', 'Uniq.',
            '-fastaout', UNIQUES_DIRECTORY + f.split('.')[0] + '.unique.fasta'
            ])

def OTUs():
    uniques = []
    for f in listdir(UNIQUES_DIRECTORY):
        if 'unique' in f:
            uniques.append(f)

    for f in uniques:
        print uniques.index(f), len(uniques)
        call(['usearch',
            '-cluster_otus', UNIQUES_DIRECTORY + f,
            '-otus', OTU_DIRECTORY + f.split('.')[0] + '.OTU.fasta',        
            '-uparseout', OTU_PARSE_DIRECTORY + f.split('.')[0] + '.up',
            '-relabel', 'OTU', '-minsize', '2'
        ])

def assign_OTUs():
    OTUs = []
    for f in listdir(OTU_DIRECTORY):
        if 'OTU' in f:
            OTUs.append(f)

    for f in OTUs:
        print OTUs.index(f), len(OTUs)
        call(['usearch',
            '-usearch_global', TRIMMED_DIRECTORY + f.split('.')[0] + '.trimmed.fasta',
            '-db', OTU_DIRECTORY + f.split('.')[0] + '.OTU.fasta',        
            '-strand', 'plus',
            '-id', '0.97',
            '-otutabout', TABLE_DIRECTORY + f.split('.')[0] + '.table.txt',
            '-biomout', TABLE_DIRECTORY + f.split('.')[0] + '.table.json',
            ''
        ])



# unzip()
# pre_join()
# L00_concat_filter()
# trim()
# unique()
# OTUs()
# assign_OTUs()
