#! /usr/bin/env python

"""
Needs modules:
    homology_builder
    beta_list
    rmsd
    rmsd_report

Use struct_align RMSD as fitness

"""


import os
from random import random, randint, choice
import random
import csv
import bisect
import collections
from operator import add
import re
from beta_list import cdf, res_choice, random_res
from itertools import imap
import operator
from homology_builder import homology_multirun, build_starter
import shutil
from prot_align import prot_align, align_rmsd, align_report, fitness_dictionary, prime_energy, prime_energy_report, prime_energy_dictionary


#############################
#####     GENERATION    #####
#############################

def config_parser():
    """
    Parse csv config file to get parameters for starting population as dict name:value
    """
    directory = os.getcwd()
    result = {}
    with open(os.path.join(directory, 'config.csv'), 'rb') as config:
        reader = csv.DictReader(config, fieldnames=['key', 'value'], delimiter='=')
        for row in reader:
            result.update({row['key']: row['value']})

    return result


def translator(res):
    d = {'A': 1, 'C': 2, 'E': 3, 'D': 4, 'G': 5,
         'F': 6, 'I': 7, 'H': 8, 'K': 9, 'M': 10,
         'L': 11, 'N': 12, 'Q': 13, 'P': 14, 'S': 15,
         'R': 16, 'T': 17, 'W': 18, 'V': 19, 'Y': 20}
    d_invert = {v:k for k,v in d.items()}

    try:
        if res in d.keys():
            res = d[res]
        elif res in d_invert.keys():
            res = d_invert[res]
    except:
        print "Wrong residue!"
    return res


def vector(base_seq, mut_factor, CDR_list):
    """base_seq - sequense in str format,
    mut_factor in range (0,1)
    return vector of integers"""
    vector = []
    list_ref = [translator(i) for i in base_seq]
    for i in range(len(list_ref)):
        if i not in CDR_list and mut_factor > random.random() and list_ref[i] not in ['C', 2]:
            vector.append(translator(res_choice()))
        else:
            vector.append(list_ref[i])
    return vector


def population_vectors(p_count, base_seq, mut_factor, CDR_list):
    """CDR insert into mutated in beta-sheet style vector"""
    population = [vector(base_seq, mut_factor, CDR_list) for i in range(p_count)]
    return population


def translator_list(list):
    d = {'A': 1, 'C': 2, 'E': 3, 'D': 4, 'G': 5,
         'F': 6, 'I': 7, 'H': 8, 'K': 9, 'M': 10,
         'L': 11, 'N': 12, 'Q': 13, 'P': 14, 'S': 15,
         'R': 16, 'T': 17, 'W': 18, 'V': 19, 'Y': 20}
    d_invert = {v:k for k,v in d.items()}

    for res in list:
        try:
            if res in d.keys():
                list[list.index(res)] = d[res]
            elif res in d_invert.keys():
                list[list.index(res)] = d_invert[res]
        except:
            print "Wrong residue!"
            list[list.index(res)] = None

    return list


def pop_creator(variants, length, min, max):
    """
    population = pop_creator(...)[i]['vec']
    """
    pop_dict = {}
    for i in range(variants):
        individual = [randint(min, max) for j in range(length)]
        pop_dict["gen_%s" % i] = {'vec': individual}
        pop_dict["gen_%s" % i].update({'rmsd': 100})
    return pop_dict


#############################
#####       FASTA       #####
#############################

def fasta_maker_dir(vec):
    """
    Outwrite separate fasta files for vectors in both: int or str kind
    """
    if isinstance(vec[0], int):
        vec_fasta = translator_list(vec)
    else:
        vec_fasta = vec
    seq = ''.join(vec_fasta)
    fasta_list = os.listdir(os.getcwd())
    fasta_dirs = [i.split('_')[-1] for i in fasta_list if i != '.directory']
#    print fasta_dirs
    try:
        num = max([int(i) for i in fasta_dirs])
    except:
        num = 0
    number = num + 1
#    print number
    os.mkdir('seq_%s' % number)
    os.chdir('seq_%s' % number)
    fasta_file = 'seq_%s.fasta' % number
    f = open(fasta_file, 'w')
    f.write('> %s \n' % number)
    f.write('%s \n' % seq)
    f.close()


def fasta_maker(pop):
    """
    Outwrite separate fasta files for population vectors in both: int or str kind
    """
    # Translate vectors if need

    counter = 0
    try:
        for vec in pop:
            if isinstance(vec[0], int):
                vec = translator_list(vec)
            else:
                pass
#            print "vector %s" % vec
            counter += 1
            seq = ''.join(vec)
#            print "%s %s" % (counter, seq)
            f = open('seq_%s.fasta' % counter, 'w')
            f.write('> %s \n' % counter)
            f.write('%s \n' % seq)
            f.close()
    except:
        pass


def grade_2(fitness_list, target):
    from operator import add
    assert(len(fitness_list) != 0)
    summed = reduce(add, (float(f) for f in fitness_list), target)
    return summed / (len(fitness_list) * 1.0)


def vec_from_fasta(seq_id):
    """ seg_id looks like seq_768 """
    fasta_file = '%s.fasta' % seq_id
    with open(fasta_file, 'r') as fl:
        data = fl.read()
        seq = data.split('\n')[1][:-1]
        vec = [translator(i) for i in seq]
        fl.close()
    return vec


#############################
#####     SELECTION     #####
#############################

def hamming_str(str1, str2):
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(imap(ne, str1, str2))


def hamming(list1, list2):
    assert len(list1) == len(list2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(list1, list2))


#############################
#####     EVOLUTION     #####
#############################

def evolve(pop, target, best, part, mutate, autbriding, cross_point):

    graded_list = []
    fitness_list = []
    fitness_min = []

#    print fitness_dictionary(coef)

    for k, v in fitness_dictionary().items():
#        print "fitness dictionary: %s - %s " % (k, v)
        fitness = v
        try:
            seq_id = k
            vector = vec_from_fasta(seq_id)     # num style
            graded_list.append((fitness, vector))
            fitness_list.append(fitness)
        except:
            pass
    print "%s scorings count" % len(fitness_list)      # Watch up scorings in pop
    print fitness_list      # Watch up scorings in pop
#    print "%s graded list: " % graded_list

    # Control the fitness average of population printing report as list and building diagramm
    fitness_report.append(grade_2(fitness_list, target))
    fitness_min.append(min(fitness_list))
    print '############## \n FITNESS MINIMUM: %s \n ################' % min(fitness_list)
    with open('fitness.txt', 'w') as config:
        config.write('FITNESS MINIMUM: %s \n ' % min(fitness_list))
        config.close()
    for fit in fitness_report:
        print "Generation %s  " % fitness_report.index(fit) + "|" * (int(float(fit) * 20) - 30)


    # take % of best scoring candidats
    graded = [x[1] for x in sorted(graded_list)]
    best_length = int(len(graded) * best)
    parents = graded[:best_length]

    # randomly add other [part] % of pop individuals to promote genetic diversity
    parent_count = int(len(pop) * part)
    while len(parents) < parent_count:
        parents.append(random.choice(graded[best_length:]))
#    print "parents: %s " % parents

    # Watch if vectors are different enough
    hamming_list = []
    for parent in parents:
        ham = hamming(parent, parents[0])
        hamming_list.append(ham)
#        print "Hamming dist: %s" % hamming_list

    # mutate some residues not touching CDRs

    for individual in parents:
        if mutate > random.random():
            pos_to_mutate = randint(0, len(individual) - 1)
            if pos_to_mutate not in CDR_list and individual[pos_to_mutate] not in ['C', 2]:       # Don`t touch CYS!
                individual[pos_to_mutate] = random.choice(range(1, 20))
#                print individual[pos_to_mutate]
     

    # crossover parents to create children
#    print "crossover parents to create children"
    parents_length = len(parents)
#    print "parents_length %s" % parents_length
    desired_length = p_count * (1 - best)
#    print "desired_length %s" % desired_length
    children = []
    counter = 0
    while len(children) < desired_length:
        counter += 1
        male = randint(0, parents_length-1)
        female = randint(0, parents_length-1)
        ham_dist = hamming(parents[male], parents[female])
#        print "%s: %s - %s %s" % (counter, male, female, ham_dist)
        if male != female and ham_dist > int(autbriding):
            male = parents[male]
            female = parents[female]
            # 1-point
            point = cross_point
            child = male[:point] + female[point:]
#           # 2-point
#            point_1 = cross_point_1
#            point_2 = cross_point_2
#            child = male[:cross_point_1] + female[cross_point_1:cross_point_2] + male[cross_point_2:]

            children.append(child)


    parents = graded[:best_length]
    parents.extend(children)
#    print "new_parents %s" % len(parents)
    with open('new_generation.csv', 'w') as ng:
        for individual in parents:
            ng.write("%s \n" % individual)

    return parents


###########################
####    MAIN CYCLE     ####
##########################


# INPUT DATA

# Sequense to evolve


# Options: take from config

input_dir = config_parser()['input_dir']
directory = os.getcwd()
autbriding = int(config_parser()['autbriding'])
base_seq = config_parser()['base_seq']
ref_file = os.path.join(input_dir, config_parser()['ref_file'])
base_mod = os.path.join(input_dir, config_parser()['base_mod'])
target = int(config_parser()['target'])
mut_factor = float(config_parser()['mut_factor'])
p_count = int(config_parser()['p_count'])
best = float(config_parser()['best'])
part = float(config_parser()['part'])
mutate = float(config_parser()['mutate'])
generations = int(config_parser()['generations'])
#coef = float(config_parser()['coef'])
cross_point = int(config_parser()['cross_point'])


# BASIC RUN

#CDR_list = range(22, 35) + range(47, 58) + range(91, len(base_seq) + 1)
CDR_list = range(20, 35) + range(47, 58) + range(91, 103)

p = population_vectors(p_count, base_seq, mut_factor, CDR_list)


parents = []
fitness_report = []
fitness_min = []


counter = 0
scoring = 10
for i in xrange(generations):
#while scoring > target:
    counter += 1
#    go = raw_input("Press Enter to start population")

# Create pop dir and copy ref_file to this directory
    root_dir = os.getcwd()
#    print "root_dir: %s" % root_dir
    pop_dir = "pop_%s" % counter
    os.mkdir(pop_dir)
    os.chdir(pop_dir)
    src = os.path.join(root_dir, '%s.mae' % ref_file)
#    print "SRC: %s" % src
    dst = os.path.join(root_dir, pop_dir)
#    print "DST: %s" % dst
    shutil.copy(src, dst)

# Write fasta files
    fasta_maker(p)
    fasta_made = len([f for f in os.listdir(os.getcwd()) if os.path.splitext(f)[1] == '.fasta'])
    print "%s fasta files made" % fasta_made

# Homology building
    build_starter()

    files = [f[:-6] for f in os.listdir(os.getcwd()) if os.path.splitext(f)[1] == '.fasta']
    for f in files:
        try:
            src = os.path.join(os.getcwd(), '%s' % f, '%s_0-out.mae' % f)
            dst = os.getcwd()
            shutil.copy(src, dst)
        except:
            pass

    mae_made = len([f for f in os.listdir(os.getcwd()) if os.path.splitext(f)[1] == '.mae']) - 1
    print "%s models build" % mae_made

# Prime Energy count
    prime_energy()
    prime_energy_made = len([f for f in os.listdir(os.getcwd()) if f.startswith('prime_energy') and os.path.splitext(f)[1] == '.csv'])
    print "Prime energy was count for %s files" % prime_energy_made
    prime_energy_report()

# Structure protein align job
    prot_align('4LLU_chainB')

# RMSD count
    align_rmsd()
    rmsd_made = len([f for f in os.listdir(os.getcwd()) if f.startswith('rmsd') and os.path.splitext(f)[1] == '.csv'])
#    print "RMSD was count for %s files" % rmsd_made
    align_report()


# Evolution run
    p = evolve(p, target, best, part, mutate, autbriding, cross_point)
    parents.append(p)
    os.chdir(root_dir)

# Write report file

    with open('fitness.txt', 'w') as config:
        config.write("pop_%s \n" % counter)
        config.write("%s models build \n" % mae_made)
        config.write("Prime energy was count for %s files \n" % prime_energy_made)
        config.close()



print "FITNESS: %s" % fitness_report

with open('fitness.txt', 'w') as config:
    config.write("FITNESS REPORT: %s" % fitness_report)
    config.close()



# Write out fitness-report




"""
# Run from previous step
#parent_pop = []
#for i in range(1, 100):
#    parent_pop.append(vec_from_fasta("seq_%s" % str(i)))
#p = parent_pop


# Save config

with open('config.txt', 'w') as config:
    config.write(info)
    config.close()

# aHER3 VHHBCD09001
base_seq_1 = "EVQLVQSGGGLVQPGGSLRLSCAASGRTSSKYAMGWFRQAPGKGTEFVATISWSDGSTYYADSVEGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAAVDVLAGTFEYEYDYWGQG"
# aHER3 VHHBCD090304
base_seq_2 = "QVQLVQSGGGLVQAGGSLRLSCAFSGRTFSMYTMGWFRQAPGKEREFVAANRGRGLSPDIADSVNGRFTISRDNAKNTLYLQMDSLKPEDTAVYYCAADLQYGSSWPQRSSAEYDYWGQGTTVTVSS"
#list_1 = [translator(i) for i in base_seq]


for pop in p:
    hamming_list = []
    ham = hamming(pop, p[0])
    hamming_list.append(ham)
    print "Hamming dist: %s" % hamming_list
"""


