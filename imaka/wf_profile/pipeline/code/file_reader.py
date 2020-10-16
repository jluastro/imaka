## file_reader.py
## Eden McEwen
## Created: October 5th 2020
## Edited: October 5th 2020
## UH REU
# This file has functions that handle reading from or searching for files

import os
import sys
import fnmatch
from pipeline.code.log_ex import *

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
    except Exception as e:
        logging.error("Exception in read_file: %s", e)
        return []

def read_file_dic(in_file, split="="):
    """ 
    Generate entries from a file
        input: in_file (readable txt file), split character default "="
        output: entries (dictonary entry per line, seperated split)
    """
    try:
        d = {}
        with open(in_file) as f:
            for line in f:
                if line[0] == "#":
                    continue
                try:
                    (key, val) = line.split(" = ")
                    d[key] = val.replace("\n", "").replace(" ", "")
                except:
                    pass
        return d
    except Exception as e:
        logging.error("Exception in read_file: %s", e)
        return {}
    
def check_target(name, target_path):
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

def read_d(entries, data_path, out_path, target_path):
    """ 
    Knowing each input is a desired directory, finds files associated
        input: elms (list of strings)
        output: d_files (list of full path to data)
                o_dirs (list of out directories for each file)
                t_files (list of target file paths for each file)
    """
    names=[]
    d_files=[]
    o_dirs=[]
    t_files=[]
    target = ""
    for d in entries:
        # check to see if this element is a new target
        t_prop = check_target(d, target_path)
        if t_prop:
            target = t_prop
        else:
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
                # create names:
                filenames = [os.path.basename(f) for f in data]
                split_out = [os.path.splitext(filename) for filename in filenames]
                name_tmp = [d+"_"+name for name, ext in split_out]
                names.extend(name_tmp)
                ## create out dir
                out_dir = out_path + d +'/'
                try:
                    os.makedirs(out_dir)
                except FileExistsError as exc:
                    pass
                o_dirs.extend([out_dir for i in data])
                t_files.extend([target for i in data])
            else:
                logging.warning("Not a dir: %s", p)
    return names, d_files, o_dirs, t_files

def read_f(entries, data_path, out_path, target_path):
    """ 
    Knowing each input is a desired files, finds files associated
        input: elms (list of strings, file paths or target names)
        output: d_files (list of full path to data)
                o_dirs (list of out directories for each file)
                t_files (list of target file paths for each file)
    """
    names=[]
    d_files=[]
    o_dirs=[]
    t_files=[]
    target = ""
    for f in entries:
        # check to see if this element is a new target
        t_prop = check_target(f, target_path)
        if t_prop:
            target = t_prop
            continue
        #does it exist?
        if os.path.isfile(f):
            ## add file to files
            d_files.append(f)
            ## create name
            filename = os.path.basename(f)
            (name, ext) = os.path.splitext(filename)
            names.append(name)
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
            logging.warning("Not a file: %s", f)
    return names, d_files, o_dirs, t_files 

    
    # !!! UNUSED !!!!
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