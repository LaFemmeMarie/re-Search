#! /usr/bin/env python

import csv
import os


def report(start):
    """
    Makes csv report from rmsd output in path directory
    """
    directory = os.getcwd()
    csv_inp = [os.path.splitext(f)[0] for f in os.listdir(directory)
               if os.path.splitext(f)[1] == '.csv' and f.startswith(start)]
    f = open('%s_report.csv' % start, 'w')
    f.write("SEQ ID, RMSD \n")
    if csv_inp == []:
        f.write("%s, %s \n" % ('No data!', 'No data!'))
    else:
        for file in csv_inp:
            reader = csv.DictReader(open(os.path.join(directory, '%s.csv' % file)))
            result = {}
            for row in reader:
                for column, value in row.iteritems():
                    result.setdefault(column, []).append(value)
            try:
                rmsd_data = result['RMSD'][0]
            except:
                rmsd_data = 'No data!'

            f.write("%s, %s \n" % (file[len(start)+1:], rmsd_data))    # cut 'start' word from seq name

    f.close()

#    return rmsd_data



def single_report(start):
    """
    Makes csv report from rmsd output in path directory
    """
    directory = os.getcwd()
    csv_inp = '%s_report.csv' % start
    reader = csv.DictReader(open(csv_inp))
    dict_start = {}
    for row in reader:
        if row['SEQ ID'] != 'No data!':
            try:
                rmsd = float((row[' RMSD ']))
                dict_start[row['SEQ ID']] = rmsd
#                print "RMSD = %s" % rmsd
            except ValueError:
#                dict_start[row['SEQ ID']] = 100
#                print "RMSD = %s" % dict_start[row['SEQ ID']]
                pass
        else:
            pass

    return dict_start


def fitness_dictionary(part):
    """
    Calculate scoring between loop and non-loop rmsd
    """
    # Try if part value is correct
    try:
        assert(part != 0)
    except(AssertionError):
        part = int(raw_input('Print part value (1-10): '))

    dict_rmsd = single_report('rmsd')
    dict_loop = single_report('loop')

    dict_report = {}
    # Difference between min and max value for rmsd in population
    dif = (max([float(v) for v in dict_rmsd.values()]) - min([float(v) for v in dict_rmsd.values()])) / part
    # add penalty as (normalised to 1 loop_rmsd) * difference/part
    for k, v in dict_rmsd.items():
#        print "%s %s" % (k,v)
#        penalty = dif * float(dict_loop[k]) / max([float(i) for i in dict_loop.values()])
        penalty = 0
        dict_report[k] = float(v) + penalty

    return dict_report


#report('loop')
#report('rmsd')

#print fitness_dictionary(4)

#single_report('loop')



