#! /usr/bin/env python

__author__ = 'lomovskayami'


import os
import schrodinger.job.jobcontrol as jc
import multiprocessing as mp
import shutil


def homology_multirun(start, end):
    """
    Homology building from separate fasta files in input from start to end file
    """
    fasta_path = os.getcwd()
    current = os.getcwd()
    runner = [os.path.join(os.environ['SCHRODINGER'], 'run'), 'homology.py']
    counter = 0
    inp = [os.path.splitext(f)[0] for f in os.listdir(fasta_path) if os.path.splitext(f)[1] == '.fasta']
    out = [f[:-10] for f in os.listdir(fasta_path) if f[-10:] == '_0-out.mae']
    input = [f for f in inp if f not in out]
#    print "run job for %s files" % len(range(start, end))

    for file in input[start:end]:
        counter += 1
        os.mkdir('%s' % file)
        os.chdir('%s' % file)
#        print('Building %s (%s of %s)' % (file, counter, len(input[start:end])))
        cmd = ['-i', os.path.join(fasta_path, '%s.fasta' % file),   # input fasta
               '-n', '1',    # Number of templates to use
               '-r',     # Re-align sequences before building the model
               '-J', '%s' % file,
#               '-b', 'sort_key=gaps',
               '-p', 'BUILD_DELETIONS=false',
               'BUILD_TRANSITIONS=false',
               'KNOWLEDGE_BASED=false'
              ]
        try:
            job = jc.launch_job(runner + cmd)
        except:
            pass

        os.chdir(current)
        src = os.path.join(current, '%s' % file, '%s_0-out.mae' % file)
        dst = current
        try:
            shutil.copy(src, dst)
        except:
            pass


def build_starter():
    """Run homology building on multirun way"""

    current = os.getcwd()
    cores = mp.cpu_count()
    directory = os.getcwd()
    inp = [os.path.splitext(f)[0] for f in os.listdir(directory) if os.path.splitext(f)[1] == '.fasta']
    length = len(inp)
    step = int(length / cores) + 1
    steps = [(step * j, step * (j + 1)) for j in range(int(cores))]


    print "Number of fasta files to build: %s" % length

    processes = [mp.Process(target=homology_multirun, args=(start, end)) for start, end in steps]
#    print "%s Homology building jobs will be run on %s cores" % (length, cores)


    for p in processes:
        p.start()

    for p in processes:
        p.join()
        os.chdir(current)



#build_starter()

#homology_multirun(0, 1)


