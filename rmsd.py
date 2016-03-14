#!/usr/bin/python

import os
import schrodinger.job.jobcontrol as jc


def rmsd_launcher(input, ref_file):
    """
    Calculate RMSD
    """
    runner = [os.path.join(os.environ['SCHRODINGER'], 'run'), 'rmsd.py']
    counter = 0
    for file in input:
        counter += 1
        print('%s %s' % (counter, file))
        cmd = ["-use_neutral_scaffold", "-c", "rmsd_%s_vs_%s.csv" % (file, os.path.basename(ref_file)), "-m",
               "-norenumber",
               "-m",
               "-verbose",
#               '-asl', "(chain.name H) AND NOT atom.ele H",
               os.path.join(INPUT_DIR, '%s_0-out.mae' % file),
               '%s.mae' % ref_file,
               '-HOST', "localhost:8"]
#        job = jc.launch_job(runner + cmd)
        try:
            job = jc.launch_job(runner + cmd)
            job.wait()
            print(job.getDuration())
            print('Finish - %s' % job.Status)
        except:
            print('Something is wrong')


def rmsd_launch(ref_file, start):
    """
    Calculate RMSD
    ref_file = "path/to/file" no ext!
    start calls the goal ASL
    """
    mae_path = os.getcwd()
    runner = [os.path.join(os.environ['SCHRODINGER'], 'run'), 'rmsd.py']

    asl_dict = {'loop': "(res.num 99-112) AND NOT atom.ele H",
                'rmsd': "not (res.num 26-32,52-57,99-124) AND NOT atom.ele H"}




    counter = 0
    for file in [f[:-10] for f in os.listdir(mae_path) if f[-10:] == '_0-out.mae']:
        counter += 1
#        print('%s %s' % (counter, file))
        cmd = ["-use_neutral_scaffold", "-c", "%s_%s.csv" % (start, file),
               "-norenumber",
               "-m",
               "-verbose",
               '-asl', "%s" % asl_dict[start],
               '%s.mae' % ref_file,
               os.path.join(mae_path, '%s_0-out.mae' % file)
               ]

#        job = jc.launch_job(runner + cmd)
        try:
            job = jc.launch_job(runner + cmd)
            job.wait()
            print(job.getDuration())
            print('Finish - %s' % job.Status)
        except:
            pass
#            print('Something is wrong')


def rmsd_individual(mae_file, ref_file):
    """
    Calculate RMSD
    """
    runner = [os.path.join(os.environ['SCHRODINGER'], 'run'), 'rmsd.py']
    cmd = ["-use_neutral_scaffold", "-c", "rmsd_%s_vs_%s.csv" % (mae_file, os.path.basename(ref_file)),
           "-norenumber",
           "-m",
           "-verbose",
#               '-asl', "(chain.name H) AND NOT atom.ele H",
           mae_file,
           ref_file,
           '-HOST', "localhost:8"]
    try:
        job = jc.launch_job(runner + cmd)
        job.wait()
        print(job.getDuration())
        print('Finish - %s' % job.Status)
    except:
        print('Something is wrong')






#rmsd_launch()

#rmsd_launcher(inp, reference_file)
#rmsd(path, reference_file)


#-use_neutral_scaffold -c test_rmsd.csv -m -norenumber -verbose
# /opt/schrodinger2015-2/run rmsd.py -use_neutral_scaffold -c test_rmsd.csv -m -norenumber -verbose vhh1_base1_mod_H3_1.mae vhh1_base1_mod_H3_4.mae


#rmsd_launch('/home/lomovskayami/schrodres/vhh/GENALG/input/vhh2_base_model_3kdm', 'loop')

#rmsd_launch('/home/lomovskayami/schrodres/vhh/GENALG/input/4LLU_chainB', 'rmsd')



