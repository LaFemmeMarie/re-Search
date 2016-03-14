#!/usr/bin/python

import os
import schrodinger.job.jobcontrol as jc
import csv
import shutil


def prot_align_all_inclusive(ref_file):
    """
    Makes prot_align job to the reference file for all mae files
    ref_file = "path/to/file" no ext!
    ref_file MUST BE IN THE CURRENT DIRECTORY!!!!!!!!!!!!!!!!!!!!!!!!
    """
    directory = os.getcwd()
    runner = [os.path.join(os.environ['SCHRODINGER'], 'utilities/structalign')]

    counter = 0
    for file in [f[:-10] for f in os.listdir(directory) if f[-10:] == '_0-out.mae']:
        counter += 1
#        print('%s %s' % (counter, file))
        cmd = ['-force',    #   force alignment even if structures are not sufficiently similar
               '%s.mae' % ref_file,
               '%s_0-out.mae' % file
               ]

#        job = jc.launch_job(runner + cmd)
        try:
            job = jc.launch_job(runner + cmd)
            job.wait()
#            print(job.getDuration())
#            print('Finish - %s' % job.Status)
        except:
            pass
#            print('Something is wrong')



def prot_align(ref_file):
    """
    Makes prot_align job to the reference file for all mae files
    ref_file = "path/to/file" no ext!
    ref_file MUST BE IN THE CURRENT DIRECTORY!!!!!!!!!!!!!!!!!!!!!!!!
    """
    directory = os.getcwd()
    runner = [os.path.join(os.environ['SCHRODINGER'], 'utilities/structalign')]

    counter = 0
    for k, v in prime_energy_dictionary().items():
        if float(v)<0:
            counter += 1
#            print('%s %s' % (counter, k))
#            print '%s_0-out.mae' % k
            cmd = ['-force',    #   force alignment even if structures are not sufficiently similar
                   '%s.mae' % ref_file,
                   '%s_0-out.mae' % k
                   ]

#            job = jc.launch_job(runner + cmd)
            try:
                job = jc.launch_job(runner + cmd)
                job.wait()
    #            print(job.getDuration())
    #            print('Finish - %s' % job.Status)
            except:
                pass
    #            print('Something is wrong')
        else:
            pass



def align_rmsd():
    """
    Takes rRMSD value rom prot_align 'rot-...' files
    """
    directory = os.getcwd()
    runner = [os.path.join(os.environ['SCHRODINGER'], 'utilities/proplister')]
    counter = 0

    for file in [f for f in os.listdir(directory) if f.startswith('rot') and f.endswith('mae')]:
        counter += 1
#        print('%s %s' % (counter, file))
        cmd = ['-p', "r_psp_StructAlign_RMSD",
               '-c', '-noheader',
               '-o', os.path.join(directory, '%s.csv' % file),
               os.path.join(directory, file)]
#        job = jc.launch_job(runner + cmd)
        try:
            job = jc.launch_job(runner + cmd)
        except:
            pass


def align_report():
    """
    Makes csv report from align_rmsd output in current directory
    """
    directory = os.getcwd()
    csv_inp = [os.path.splitext(f)[0] for f in os.listdir(directory)
               if os.path.splitext(f)[1] == '.csv' and f.startswith('rot')]
    f = open('align_report.csv', 'w')
    f.write("SEQ ID,RMSD\n")
    if csv_inp == []:
        f.write("%s,%s\n" % ('No data!', 'No data!'))
    else:
        for file in csv_inp:
            reader = csv.reader(open(os.path.join(directory, '%s.csv' % file)))
            for row in reader:
                try:
                    rmsd_data = float(row[0])
                except:
                    rmsd_data = 'No data!'

                f.write("%s,%s\n" % (file.split('_0-out')[0][4:], rmsd_data))    # cut 'rot-' word from seq name
    f.close()
#    return rmsd_data


def prime_energy():
    """
    Takes PrimeEnergy value from homology_building files ENERGY_BASED mode
    """
    directory = os.getcwd()
    runner = [os.path.join(os.environ['SCHRODINGER'], 'utilities/proplister')]
    counter = 0

    for file in [f for f in os.listdir(directory) if f.endswith('_0-out.mae')]:
        counter += 1
#        print('%s %s' % (counter, file))
        cmd = ['-p', "r_psp_Prime_Energy",
               '-c', '-noheader',
               '-o', os.path.join(directory, 'prime_energy_%s.csv' % file.split('_0-out')[0]),
               os.path.join(directory, file)]
#        job = jc.launch_job(runner + cmd)
        try:
            job = jc.launch_job(runner + cmd)
        except:
            pass


def prime_energy_report():
    """
    Makes csv report from csv output from proplister prime_energy in current directory
    """
    directory = os.getcwd()
    csv_inp = [os.path.splitext(f)[0] for f in os.listdir(directory)
               if os.path.splitext(f)[1] == '.csv' and f.startswith('prime_energy')]
    f = open('prime_energy_report.csv', 'w')
    f.write("SEQ ID,prime_energy\n")
    if csv_inp == []:
        f.write("%s,%s\n" % ('No data!', 'No data!'))
    else:
        for file in csv_inp:
            reader = csv.reader(open(os.path.join(directory, '%s.csv' % file)))
            for row in reader:
                try:
                    prime_energy = float(row[0])
                except:
                    prime_energy = 'No data!'

                f.write("%s,%s\n" % (file.split('energy_')[-1], prime_energy))    # cut 'prime_energy' word from seq name
    f.close()
#    return rmsd_data


def prime_energy_dictionary():
    """
    Calculate scoring from rmsd report
    """
    directory = os.getcwd()
    csv_inp = 'prime_energy_report.csv'
    reader = csv.DictReader(open(csv_inp))
    prime_energy_dict = {}
    for row in reader:
        if row['SEQ ID'] != 'No data!':
            try:
                prime_energy = float((row['prime_energy']))
                prime_energy_dict[row['SEQ ID']] = prime_energy
#                print "Prime Energy = %s" % prime_energy
            except ValueError:
                pass
        else:
            pass

    return prime_energy_dict


def fitness_dictionary():
    """
    Calculate scoring from rmsd report
    """
    directory = os.getcwd()
    csv_inp = 'align_report.csv'
    reader = csv.DictReader(open(csv_inp))
    rmsd_dict = {}
    for row in reader:
        if row['SEQ ID'] != 'No data!':
            try:
                rmsd = float((row['RMSD']))
                rmsd_dict[row['SEQ ID']] = rmsd
#                print "RMSD = %s" % rmsd
            except ValueError:
                pass
        else:
            pass

    return rmsd_dict


def proplister(file, property):
    """
    Takes property value from X-type files and puts it to .csv file
    output: csv (just value)
    input: file == FULL NAME with EXT!!!
    """
    runner = [os.path.join(os.environ['SCHRODINGER'], 'utilities/proplister')]
    counter = 0

    cmd = ['-p', "%s" % property,
           '-c', '-noheader',
           file]
    job = jc.launch_job(runner + cmd)
    job.wait()




#prot_align('4LLU_chainB')
#align_rmsd()
#align_report()

#prime_energy()
#prime_energy_report()
#print prime_energy_dictionary()
#print fitness_dictionary()

#prot_align('4LLU_chainB.mae')


# Console commands:
# /opt/schrodinger2015-2/utilities/structalign goal_file.mae file1.mae
#  $SCHRODINGER/utilities/proplister -p 'r_psp_StructAlign_RMSD' -c -noheader rot-file.mae
#  $SCHRODINGER/utilities/proplister -noheader -c -p 'r_psp_Prime_Energy' lambda22/pop_1/seq_14_0-out.mae




