#!/usr/bin/python

import os
import csv


def fitness_list():
    """
    Makes csv report from rmsd output in path directory
    """
    directory = os.getcwd()

    csv_inp = 'rmsd_report.csv'
    reader = csv.DictReader(open(csv_inp))
    fitness = []
    for row in reader:
        if row['SEQ ID'] != 'No data!':
            try:
                rmsd = float((row[' RMSD ']))
            except:
                pass
            finally:
                fitness.append(rmsd)
        else:
            pass

    return fitness


def grade(fitness_list, target=0):
    from operator import add
    assert(len(fitness_list) != 0)
    fit = [f for f in fitness_list if f != "No data!"]
    summed = reduce(add, (float(f) for f in fit), target)
    return summed / (len(fit) * 1.0)

print grade(fitness_list())



"""
def fitness_pict():

    directory = os.getcwd()
#    pop_list = [i.split('_')[-1] for i in os.listdir(directory) if i.startswith('pop_')]
    pop_list = [i for i in os.listdir(directory) if i.startswith('pop')]
    for i in pop_list:
        csv_inp = os.path.join(i, 'rmsd_report.csv')
        reader = csv.DictReader(open(csv_inp))
        fitness = []
        for row in reader:
            if row['SEQ ID'] != 'No data!':
                try:
                    rmsd = float((row[' RMSD ']))
                except:
                    pass
                finally:
                    fitness.append(rmsd)
            else:
                pass

    return fitness




for i in fitness_report:
        print "Generation %s  " % fitness_report.index(i) + "|" * int(i)

"""