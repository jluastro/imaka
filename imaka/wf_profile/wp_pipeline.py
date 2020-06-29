## Eden McEwen
## Created: June 15th 2020
## Edited: June 25th 2020
## UH REU
# This file takes in a text file of dates, finds them on ehu, and then returns correlation outputs
import os
import sys
import fnmatch

from data_pipeline import *


def data_proc(name, data_f, param_f, target_f, out_d):
    # create an object
    print("======== %s ========" % name)
    curr_data = DataPipe(name, data_f, out_d, tmax=250)
    curr_data.set_target(target_f)
    try:
        # generate acor
        print("===> acor")
        t0 = time.time()
        curr_data.acor_gen()
        t1 = time.time()
        print("<=== acor done, %s s"% str(t1-t0))
        # generate ccor
        print("===> ccor")
        t2 = time.time()
        curr_data.ccor_gen()
        t3 = time.time()
        print("<=== ccor done, %s s"% str(t3-t2))
        # create fits
        print("===> writing fits")
        out = curr_data.fits_write()
        print("<=== %s"% out)
        print("===> Graphing")
        # graph acor
        curr_data.acor_graph(t_list=[0,5,10,20,30], avg_sub=True, avg_len=5)
        curr_data.acor_animate_avg(dt_max=200, avg_sub=True, avg_len=5)
        # graph ccor
        curr_data.ccor_graph_all(avg_sub=True, avg_len=5)
        curr_data.cor_animate_all(10, avg_sub=True, avg_len=5)
        print("===> complete")
    except:
        print("ran into error")

def param_dict(param_f):
    #TODO
    if not os.path.isfile(param_f):
        return None
    f = open(param_f, "r")
    if f.mode == 'r':
        contents = f.read()
    return
        
def format_e(elms):
    """ 
    Takes in a list of elements, converts to float or int if possible
    input: elms=list of strings
    output: elms=list of strings, ints, floats
    """
    for i, e in enumerate(elms):
        try:
            val = int(e)
            elms[i] = val
        except ValueError:
            try:
                val = float(e)
                elms[i] = val
            except ValueError:
                elms[i] = e
    return elms if len(elms)>1 else elms[0]

def read_file(in_file):
    """ 
    Generate entries from a file
    input: in_file=readable txt file
    output: entries=list of strings, seperated by line
    """
    try:
        with open(in_file,'r') as i:
            entries = i.readlines() 
            entries = [e.replace('\n', '') for e in entries]
        return entries
    except (OSError, IOError) as e:
        print("Error: opening input file")
        return []


if __name__ == '__main__':
    in_type = sys.argv[1]
    in_file = sys.argv[2]
    FILE = False
    data_path = "/home/imaka/data/"
    out_path = "/home/emcewen/out/"
    d_files=[]
    params=[]
    o_dirs=[]

    ##### INPUT 1: type of txt file #####
    # for each line check if it is a valid directory or file 
    if in_type=='-d'or in_type=='-dir':
        FILE = False
    elif in_type=='-f'or in_type=='-file':
        FILE = True
    else:
        print("Error: file type not recognized")
        sys.exit(0)

    ##### INPUT 2: entry txt file #####
    #look at infile, error if infile 
    entries = read_file(in_file)
    
    ### Building list of inputs, output directories
    if FILE:
        # look at each entry to see if exitst
        for e in entries:
            #does it exist?
            if os.path.isfile(e):
                ## add file to files
                d_files.append(e)
                ## find imaka params
                p = os.path.dirname(e)
                parm = p + '/imakaparm.txt'
                if os.path.isfile(parm):
                    params.append(parm)
                else:
                    params.append(False)
                ## create out dir
                if data_path in p: p.replace(data_path, '')
                if 'ao/' in p: p.replace('ao/', '')
                out_dir = out_path + p
                try:
                    os.mkdir(out_dir)
                except FileExistsError as exc:
                    pass
                o_dirs.append(out_dir)
            else:
                print("Not a file: " + e)
    else:
        #create a list of all the files usable
        for d in entries:
            p = data_path + d + '/ao/'
            if os.path.isdir(p):
                # list of files
                files = os.listdir(p)
                ## find and add all data files
                data =  [p + fn for fn in files if fnmatch.fnmatch(fn, 'ao*o.fits')]
                print('data: ', data)
                d_files.extend(data)
                ## find imaka params
                parm = p + '/imakaparm.txt'
                if os.path.isfile(parm):
                    params.extend([parm for i in data])
                else:
                    params.extend([False for i in data])
                ## create out dir
                out_dir = out_path + d +'/'
                try:
                    os.makedirs(out_dir)
                except FileExistsError as exc:
                    pass
                o_dirs.extend([out_dir for i in data])
            else:
                print("Not a dir: " + p)
        #find imaka params
    print('files: ', d_files)
    target= "/home/emcewen/data/target_input/Orion2_USNOB_GS.txt"
    #### Applying Pipeline code ####
    for i, f in enumerate(d_files):
        print("============ Starting file %s of %s"%(i,len(d_files)))
        filename = os.path.basename(f)
        (name, ext) = os.path.splitext(filename)
        data_f = f
        param_f = params[i]
        out_d = o_dirs[i]
        print("name: %s"% name)
        print("d_file: %s"% data_f)
        print("o_dirs: %s"% out_d)
        data_proc(name, data_f, param_f, target, out_d)
        
    