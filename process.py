from os import listdir
import os
import sys
from subprocess import call
from os.path import expanduser
import json

ZIPPED_DIRECTORY = "0_zipped/"
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
MIN_LENGTH_THRESHOLD = 140

FULL_FASTA = 'full'

walk_dir = ''
print sys.argv
if len(sys.argv) == 1:
    walk_dir = ZIPPED_DIRECTORY
else:
    walk_dir = sys.argv[1]

def makeFolders(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
makeFolders(folders)

prefix = 'ssr_'

def unzip(lookup_table):
    for root, subdirs, files in os.walk(walk_dir):
        if len(files) > 0:
            for f in files:
                if f.split('.')[-1] == 'gz':
                    identifier = f[9:]
                    kID = os.path.join(root, f).split('/')[-2]
                    if kID not in lookup_table:
                        print 'kID ' + str(kID) + ' not found!!!'
                        break
                    print 'Copying/Unzipping ' + f                   
                    call(["cp", os.path.join(root, f), UNZIPPED_DIRECTORY])
                    newname = lookup_table[kID] + identifier
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


        command = 'usearch -fastq_filter ' + CONCAT_DIRECTORY + pre + 'L001.concat.fastq -relabel @ -fastq_minlen ' + str(MIN_LENGTH_THRESHOLD) + ' -fastqout ' + CONCAT_DIRECTORY + pre + 'L001.filtered.fastq'
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
        call(['vsearch',
            '--threads', '16',
            '--derep_fulllength', TRIMMED_DIRECTORY + f,
            '--sizeout',
            '--output', UNIQUES_DIRECTORY + f.split('.')[0] + '.unique.fasta',
            '--relabel', f.split('_')[0] + '.',
            '--fasta_width', '0'
            ])

def full_concatentation_all_reads():
    trimmed = []
    for f in listdir(TRIMMED_DIRECTORY):
        if 'trimmed' in f:
            trimmed.append(TRIMMED_DIRECTORY + f)

    command = 'cat '
    for c in trimmed:
        command += c
        command += ' '
    command += '> all_reads.fasta'
    # print 'concatenating', trimmed
    call(command, shell=True)
    # print command

def full_concatentation_unique():
    uniques = []
    for f in listdir(UNIQUES_DIRECTORY):
        if 'unique' in f and FULL_FASTA not in f:
            uniques.append(UNIQUES_DIRECTORY + f)

    command = 'cat '
    for c in uniques:
        command += c
        command += ' '
    command += '> ' + FULL_FASTA + '.unique.fasta'
    # print 'concatenating', trimmed
    call(command, shell=True)
    # print command

def dereplicate_full_uniques():
    call(['vsearch',
        '--threads', '16',
        '--derep_fulllength', 'full.unique.fasta',
        '--minuniquesize', '2',
        '--sizein',
        '--sizeout',
        '--fasta_width', '0',
        '--uc', 'full.derep.uc',
        '--output', 'full.derep.fasta'
        ])

def precluster():
    call(['vsearch',
        '--threads', '16',
        '--cluster_size', 'full.derep.fasta',
        '--id', '.98',
        '--strand', 'plus',
        '--sizein',
        '--sizeout',
        '--fasta_width', '0',
        '--uc', 'full.preclustered.uc',
        '--centroids', 'full.preclustered.fasta'
        ])

def denovo_chimera():
    call(['vsearch',
        '--threads', '16',
        '--uchime_denovo', 'full.preclustered.fasta',
        '--sizein',
        '--sizeout',
        '--fasta_width', '0',
        '--nonchimeras', 'full.denovo.nonchimeras.fasta'
        ])

def reference_chimera():
    call(['vsearch',
        '--threads', '16',
        '--uchime_ref', 'full.denovo.nonchimeras.fasta',
        '-db', 'silva.bacteria/silva.gold.ng.fasta',
        '--sizein',
        '--sizeout',
        '--fasta_width', '0',
        '--nonchimeras', 'full.ref.nonchimeras.fasta'
        ])

def extract_nonchimeras():
    call('perl map.pl full.derep.fasta full.preclustered.uc full.ref.nonchimeras.fasta > full.nonchimeras.derep.fasta', shell = True)
    call('perl map.pl full.unique.fasta full.derep.uc full.nonchimeras.derep.fasta > full.nonchimeras.fasta', shell = True)

def OTUs():
    call(['vsearch',
        '--threads', '16',
        '--cluster_size', 'full.nonchimeras.fasta',
        '--id',' 0.97',
        '--strand', 'plus',
        '--sizein',
        '--sizeout',
        '--fasta_width', '0',
        '--uc', 'full.clustered.uc',
        '--relabel', 'OTU_',
        '--centroids', 'full.otus.fasta',
        '--otutabout', 'full.otutab.txt'
        ])
 
# lookup_table = json.loads(open('lookup_table.json','r').read())


# unzip(lookup_table)
# pre_join()
# L00_concat_filter()
# trim()
# unique()
# full_concatentation_unique()
# dereplicate_full_uniques()
# precluster()
# denovo_chimera()
# reference_chimera()
# extract_nonchimeras()
# OTUs()
# full_concatentation_all_reads()
