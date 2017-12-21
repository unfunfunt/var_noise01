#!/usr/bin/env python

import spider_analysis as sa
import os
import sys
import argparse
import subprocess
from spider_analysis.batch import qsub_group, qsub


prod = 'var_noise01'
prod_dir = 'varnoise01dir'

parser = argparse.ArgumentParser()
parser.add_argument('--fpu', action='store', type=int, choices=range(1,7), help='The fpu to use')
parser.add_argument('--column', action='store', type=int, choices=range(0, 16), help='the column to use')
parser.add_argument('--row', action='store', type=int, choices=range(0, 33), help='the row to use')
parser.add_argument('--cpu-speed', type=float, default=1.0, help='relative cpu speed factor')
parser.add_argument('--queue', help='Queue to submit the jobs to')
parser.add_argument('--nodes', type=str, default=1, help='number of nodes to submit the job to')
parser.add_argument('--exclude', type=str, default=None, nargs='+', help='nodes to exclude from the job')
parser.add_argument('--env-script', help='environment setup script')
parser.add_argument('--slurm', action='store_true', default=False, help='use SLURM scheduling instead of PBS')
parser.add_argument('--group', type=int, default=1, help='Number of processes per job')
parser.add_argument('--out-path', default=prod, help='Output directory')
parser.add_argument('--use-cput', action='store_true', help='use cpu time')

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--cluster', action='store_true', help='run on a cluster')
group.add_argument('--bee', action='store_true', default=False, help='Run on bee')
group.add_argument('--slow', action='store_true', default=False, help='Run the code serially')
group.add_argument('--subprocess', action='store_true', default=False, help="Don't set this from the command line")

args = parser.parse_args()

#product name
product = 'var_noise01'

testing = True

data_root = os.environ['SPIDER_DATA_ROOT']
if os.path.isdir(data_root) is False:
    sys.exit("SPIDER_DATA_ROOT not set or isn't a directory")

U = sa.Unifile(default_latest=True, verbose=3)

if not args.subprocess:
    opts = {'unifile_version' : 5,
            'ruleset_version' : 3,
            'creator' : 'MNG',
            'deps' : ['dcclean08', 'stepstitch07', 'planck02', 'rwn01', 'yssn17', 'point06'],
            'aux_deps' : ['microchunks05'],
            'testing' : testing,
            'encoding' : 'sie',
            'date' : (2017, 10, 18)}

    P = U.new_product(product, **opts)

    compiler = os.environ.get('CC', 'gcc')
    cflags = os.environ.get('CFLAGS', '-O3') + ' -std=c99 '

    source_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), prod_dir)

    compile_cmd = compiler + ' ' + cflags + ' -o ' + prod + ' ' + os.path.join(source_dir, 'noisepowspec.c ') + os.path.join(source_dir, 'psdfftw3.c ') + os.path.join(source_dir, 'gslfits.c ') + os.path.join(source_dir, 'lsperiodigram.c ') + os.path.join(source_dir, 'dipole.c') + ' -lgetdata -lgsl -lqpoint -lchealpix_qp -fopenmp -lsofa_c -lslarefro -lgslcblas -lm -lfftw3 -l:libnfft3.so'

    print 'compiling: ' + compile_cmd
    success = os.system(compile_cmd)
    if success is not 0:
        sys.exit("failed to compile", product)


    #make the dirfile with all the fields
    if testing is False:
        for chan in U.get_channel_set():
            for fieldname in ['_a', '_b', '_c', '_d', '_x']:
                scaling = 10/32768.0
                if fieldname == '_x':
                    scaling = 100/32768.0
                if fieldname == '_b':
                    scaling = 0.1/32768.0

                P.add_raw_field(chan+fieldname, dtype='int16', isupper=False, spf=20, sync=False)
                P.add_lincom_field((chan + fieldname).upper(), chan+fieldname, scaling, 0, overwrite=True, sync=False, isupper=False)

    P.close()
    U.close()

    print 'done adding fields'

#run the product

if args.fpu is None:
    fpu_list = range(1, 7)
else:
    fpu_list = [args.fpu]

if args.column is None:
    column_list = range(0, 16)
else:
    column_list = [args.column]

if args.row is None:
    row_list = range(0, 33)
else:
    row_list = [args.row]

if args.slow:
    for fpu in fpu_list:
        for row in row_list:
            for column in column_list:
                command = './' + prod + ' -r -y x' + str(fpu) + ' c' + str(column) + ' r' +str(row)
                success = os.system(command)
                if success is not 0:
                    sys.exit('failed to run command ' + command)

elif args.bee:
    #do something to make it use some cores on bee
    pids = set()
    for fpu in fpu_list:
        arglist = ['./'+ prod +'.py', '--fpu', str(fpu), '--subprocess']
        if args.column is not None:
            arglist.extend(['--column', str(args.column)])
        if args.row is not None:
            arglist.extend(['--row', str(args.row)])
        p = subprocess.Popen(arglist)
        pids.add(p.pid)
    while pids:
        pid, retval=os.wait()
        print('process ' + str(pid) + ' finished')
        pids.remove(pid)

elif args.cluster:
    #write submit scripts and stuff
    script_path = os.path.join(os.getenv('SPIDER_TOOLS_PATH'), 'analysis', 'generator_scripts')
    path2log = os.path.join(script_path, '{}_logs'.format(prod))

    job_opts = {
            'nodes' : args.nodes,
            'exclude': args.exclude,
            'ppn': args.group,
            'queue': args.queue,
            'workdir': path2log,
            'omp_threads': 1,
            'env_script': args.env_script,
            'delete': False
        }

    cpu_time = 24 # hours

    if args.use_cput:
        job_opts['cput'] = cpu_time * args.group/ args.cpu_speed
    else:
        job_opts['wallt'] = cpu_time / args.cpu_speed

    if args.slurm:
        job_opts['scheduler'] = 'slurm'

    jobs_list = []
    for fpu in fpu_list:
        if args.column is None:
            for col in column_list:
                job_str = os.path.join(script_path, '{}.py'.format(prod)) + ' --fpu {} --col {} --subprocess'.format(fpu, col)
                if args.column is not None:
                    job_str += '--column {}'.format(args.column)
                jobs_list.append(job_str)

        else:
            for row in row_list:
                job_str = os.path.join(script_path, '{}.py'.format(prod)) + ' --fpu {} --row {} --subprocess'.format(fpu, row)
                jobs_list.append(job_str)

    job_opts['name'] = prod
    job_ids = qsub_group(jobs_list, args.group, **job_opts)

else:
    #run it 
    for opts in U.event_partition("lst_day", event="all"):
        start, count = U.get_frame_count(return_start=True, **opts)
        command = './' + prod + ' -r -y x' + str(args.fpu) 
        if args.column is not None:
            command += ' c' + str(args.column) 
        if args.row is not None:
            command += ' r' + str(args.row)
        command += ' -s ' + str(20*start) + ' -e ' + str(20*(start + count)) + ' >> ./'+ prod_dir +'/x' + str(args.fpu) + '.log'
        success = os.system(command)
        print command
        if success is not 0:
            sys.exit('error running ' + command)
