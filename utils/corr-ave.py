from __future__ import print_function
import re
import argparse
import numpy as np
import h5py
import sys
import time

def get_root(fnames):
    # Return root of HDF5 files with names in fnames. If the root is
    # not unique, throw an exception
    roots = []
    for fn in fnames:
        with h5py.File(fn, "r") as fp:
            roots += list(fp.keys())
    assert len(set(roots)) == 1, "Files do not have unique root"
    return roots[0]

def get_src_pos(fnames, root):
    # Return a list of source positions from HDF5 files with names in
    # fnames. These should be the names of the second-level group. If
    # the group names are not as expected, throw an assertion exception
    spos = []
    for fn in fnames:
        with h5py.File(fn, "r") as fp:
            spos += list(fp[root].keys())
    assert len(set(spos)) == len(fnames), "Number of unique source positions != len(fnames)"
    s = "sx[0-9]{2}sy[0-9]{2}sz[0-9]{2}st[0-9]{2}"
    for spo in spos:
        assert re.match(s, spo) is not None, "One or more second-level groups not of the form %s" % s
    return spos

def get_dset_names(fnames, root, spos):
    # Returns the full list of dataset names in fnames, under
    # root/spos These must be unique for all files, otherwise an
    # assertion error is thrown
    names = {fn: list() for fn in fnames}
    for fn,sp in zip(fnames,spos):
        with h5py.File(fn, "r") as fp:
            grp = fp[root][sp]
            grp.visititems(lambda x,t: names[fn].append(x) if type(t) is h5py.Dataset else None)
    for fn in fnames[1:]:
        assert names[fnames[0]] == names[fn], "File structure missmatch"
    names = names[fnames[0]]
    # Strip dataset name from names. Remember we have two datasets per
    # last-level group: arr and mvec
    names = ["/".join(n.split("/")[:-1]) for n in names if 'arr' in n]
    return names

def get_dset(fname, name):
    # Returns the dataset with name "name" of fname
    with h5py.File(fname, "r") as fp:
        dset = np.array(fp[name])
    return dset

def init_h5file(fname, root, attrs=None):
    with h5py.File(fname, "w") as fp:
        fp.create_group(root)
        if attrs is not None:
            for k,v in attrs.items():
                fp[root].attrs[k] = np.array(v,'S')
    return

def write_dset(fname, grp_name, arr):
    with h5py.File(fname, "a") as fp:
        grp = fp.create_group(grp_name)
        grp.create_dataset('arr', arr['arr'].shape, dtype=arr['arr'].dtype, data=arr['arr'])
        grp.create_dataset('mvec', arr['mvec'].shape, dtype=arr['mvec'].dtype, data=arr['mvec'])
    return        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fn0", type=str, metavar="FNAME", help=argparse.SUPPRESS)
    parser.add_argument("fn1", type=str, metavar="FNAME", nargs="+", help="two or more file names")
    parser.add_argument("-o", "--output", metavar="F", type=str, default='ave.h5',
                        help="output file name (default: ave.h5)")
    args = parser.parse_args()
    fnames = [args.fn0] + args.fn1
    output = args.output    
    root = get_root(fnames)
    spos = get_src_pos(fnames, root)
    names = get_dset_names(fnames, root, spos)
    data = {sp: {} for sp in spos}
    for fn,sp in zip(fnames, spos):
        for n in names:
            ds = root + "/" + sp + "/" + n + "/arr"
            arr = get_dset(fn, ds)
            ds = root + "/" + sp + "/" + n + "/mvec"
            mvec = get_dset(fn, ds)
            data[sp][n] = {"arr": arr, "mvec": mvec}
    # Now check if the momentum vectors of each dataset match
    for n in names:
        for sp in spos[1:]:
            a = data[spos[0]][n]["mvec"].tolist()
            b = data[sp][n]["mvec"].tolist()
            assert a == b, " Missmatch in momentum vectors"
    # Average arr over source positions Momentum vectors are taken
    # from first file. We've allready checked they are identical over
    # files
    ave = {n: {'arr': np.array([data[sp][n]['arr'] for sp in spos]).mean(axis=0),
               'mvec': data[spos[0]][n]['mvec']}
           for n in names}
    nsrc = len(fnames)
    top = root + "/ave%d" % nsrc
    init_h5file(output, top, attrs={"Source positions": spos})
    for n in names:
        grp = top + "/" + n
        write_dset(output, grp, ave[n])
    return 0

if __name__ == "__main__":
    sys.exit(main())
