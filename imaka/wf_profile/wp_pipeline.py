## Eden McEwen
## Created: June 15th 2020
## Edited: June 25th 2020
## UH REU
# This file takes in a text file of dates, finds them on ehu, and then returns correlation outputs
import os
import sys
import fnmatch

from data_pipeline import *

# Global variables, change when running of diff systems
data_path = "/home/imaka/data/"
out_path = "/home/emcewen/out/"
target_path = "/home/emcewen/data/target_input/"


def data_proc(name, data_f, out_d, target_f, s_sub=False, tt_sub = False):
    # create an object
    print("======== %s ========" % name)
    curr_data = DataPipe(name, data_f, out_d, tmax=200)
    curr_data.set_target(target_f)
    curr_data.set_ssub(s_sub)
    curr_data.set_ttsub(tt_sub)
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
    return None
        
def read_file(in_file):
    """ 
    Generate entries from a file
        input: in_file (readable txt file)
        output: entries (list of strings, seperated by line)
    """
    try:
        with open(in_file,'r') as i:
            entries = i.readlines() 
            entries = [e.replace('\n', '') for e in entries]
        return entries
    except (OSError, IOError) as e:
        print("Error: opening input file")
        return []

def format_e(elms):
    """ 
    Takes in a list of elements, converts to float or int if possible
        input: elms (list of strings)
        output: elms (list of strings, ints, floats)
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

def check_target(name):
    """
    Takes in the target name, returns the full path to file if there is one, otherwise false
        input: string
        output: path (string) or False
    """
    prefix = target_path
    suffix = "_USNOB_GS.txt"
    file = prefix + name + suffix
    if os.path.isfile(file):
        return file
    else:
        return False

def read_d(entries):
    """ 
    Knowing each input is a desired directory, finds files associated
        input: elms (list of strings)
        output: d_files (list of full path to data)
                o_dirs (list of out directories for each file)
                t_files (list of target file paths for each file)
    """
    d_files=[]
    o_dirs=[]
    t_files=[]
    target = ""
    for d in entries:
        # check to see if this element is a new target
        t_prop = check_target(d)
        if t_prop:
            target = t_prop
            continue
        # if not a target file, will check to see if the dir path is available
        p = data_path + d + '/ao/'
        if os.path.isdir(p):
            # list of files
            files = os.listdir(p)
            ## find and add all data files
            data =  [p + fn for fn in files if fnmatch.fnmatch(fn, 'ao*o.fits')]
            #print('data: ', data)
            if not data:
                data =  [p + fn for fn in files if fnmatch.fnmatch(fn, 'ao*.fits')]
            d_files.extend(data)
            ## create out dir
            out_dir = out_path + d +'/'
            try:
                os.makedirs(out_dir)
            except FileExistsError as exc:
                pass
            o_dirs.extend([out_dir for i in data])
            t_files.extend([target for i in data])
        else:
            print("Not a dir: " + p)
    return d_files, o_dirs, t_files

def read_f(entries):
    """ 
    Knowing each input is a desired files, finds files associated
        input: elms (list of strings, file paths or target names)
        output: d_files (list of full path to data)
                o_dirs (list of out directories for each file)
                t_files (list of target file paths for each file)
    """
    d_files=[]
    o_dirs=[]
    t_files=[]
    target = ""
    for f in entries:
        # check to see if this element is a new target
        t_prop = check_target(f)
        if t_prop:
            target = t_prop
            continue
        #does it exist?
        if os.path.isfile(f):
            ## add file to files
            d_files.append(f)
            ## create out dir
            if data_path in p: p.replace(data_path, '')
            if 'ao/' in p: p.replace('ao/', '')
            out_dir = out_path + p
            try:
                os.mkdir(out_dir)
            except FileExistsError as exc:
                pass
            o_dirs.append(out_dir)
            t_files.append(target)
        else:
            print("Not a file: " + f)
    return d_files, o_dirs, t_files




###################################################################
##################### MAIN FUNCTION ###############################
###################################################################

if __name__ == '__main__':
    in_type = sys.argv[1]
    in_file = sys.argv[2]
    FILE = False
    d_files=[]
    o_dirs=[]
    t_files=[]
    
    ##### INPUT 2: entry txt file #####
    #look at infile, error if infile 
    entries = read_file(in_file)

    ##### INPUT 1: type of txt file #####
    # for each line check if it is a valid directory or file 
    if in_type=='-d':
        d_files, o_dirs, t_files = read_d(entries)
    elif in_type=='-f':
        d_files, o_dirs, t_files = read_f(entries)
    else:
        print("Error: file type not recognized")
        sys.exit(0)

    print('files: ', d_files)
    #print('targets: ', t_files)
    #### Applying Pipeline code ####
    for i, f in enumerate(d_files):
        print("============ Starting file %s of %s"%(i + 1,len(d_files)))
        filename = os.path.basename(f)
        (name, ext) = os.path.splitext(filename)
        data_f = f
        out_d = o_dirs[i]
        target_f = t_files[i]
        print("name: %s"% name)
        print("d_file: %s"% data_f)
        print("o_dirs: %s"% out_d)
        print("t_file: %s"% target_f)
        data_proc(name, data_f, out_d, target_f)
        data_proc(name, data_f, out_d, target_f, s_sub=True)
        data_proc(name, data_f, out_d, target_f, s_sub=True, tt_sub=True)
        

########### Old code for finding imaka params, not used
### for -d
## find imaka params
#                parm = p + '/imakaparm.txt'
#                if os.path.isfile(parm):
#                    params.extend([parm for i in data])
#                else:
#                    params.extend([False for i in data])

#### for -f 
## find imaka params
#                p = os.path.dirname(e)
#                parm = p + '/imakaparm.txt'
#                if os.path.isfile(parm):
#                    params.append(parm)
#                else:
#                    params.append(False)