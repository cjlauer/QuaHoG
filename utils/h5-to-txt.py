from __future__ import print_function
import argparse
import numpy as np
import h5py
import sys
import time
import re

ND = 4
NS = 4

meson_gams = ["=1=", "=g5=", "=gx=", "=gy=", "=gz=", "=gt=", "=g5gx=", "=g5gy=", "=g5gz=", "=g5gt=",]
meson_flvs = ["up", "dn"]
nucl_flvs = ["ppm", "pmm"]
nucl_intrps = ["1-1", "1-2", "2-1", "2-2"]
thrp_gams = ["=loc:1=", "=loc:g5=", "=loc:g0=", "=loc:gx=", "=loc:gy=", "=loc:gz=",
             "=loc:g5g0=", "=loc:g5gx=", "=loc:g5gy=", "=loc:g5gz=",
             "=loc:g5si0x=", "=loc:g5si0y=", "=loc:g5si0z=",
             "=loc:g5sixy=", "=loc:g5sixz=",
             "=loc:g5siyz=",
             "=noe:g0=", "=noe:gx=", "=noe:gy=", "=noe:gz=",
             "=der:g0D0:sym=", "=der:gxD0:sym=", "=der:gyD0:sym=", "=der:gzD0:sym=",
             "=der:gxDx:sym=", "=der:gyDx:sym=", "=der:gzDx:sym=",
             "=der:gyDy:sym=", "=der:gzDy:sym=",
             "=der:gzDz:sym=",
             "=der:g5g0D0:sym=", "=der:g5gxD0:sym=", "=der:g5gyD0:sym=", "=der:g5gzD0:sym=",
             "=der:g5gxDx:sym=", "=der:g5gyDx:sym=", "=der:g5gzDx:sym=",
             "=der:g5gyDy:sym=", "=der:g5gzDy:sym=",
             "=der:g5gzDz:sym=",
             "=der:g5gxD0:asy=", "=der:g5gyD0:asy=", "=der:g5gzD0:asy=",
             "=der:g5gyDx:asy=", "=der:g5gzDx:asy=",
             "=der:g5gzDy:asy=",
             "=der:g5si00Dx:sym=", "=der:g5si00Dy:sym=", "=der:g5si00Dz:sym=",
             "=der:g5si0xDx:sym=", "=der:g5si0xDy:sym=", "=der:g5si0xDz:sym=",
             "=der:g5si0yDy:sym=", "=der:g5si0yDz:sym=", "=der:g5si0zDz:sym=",
             "=der:g5six0D0:sym=", "=der:g5six0Dx:sym=", "=der:g5six0Dy:sym=",
             "=der:g5six0Dz:sym=", "=der:g5sixxDy:sym=", "=der:g5sixxDz:sym=",
             "=der:g5sixyDy:sym=", "=der:g5sixyDz:sym=", "=der:g5sixzDz:sym=",
             "=der:g5siy0D0:sym=", "=der:g5siy0Dx:sym=", "=der:g5siy0Dy:sym=",
             "=der:g5siy0Dz:sym=", "=der:g5siyxDx:sym=", "=der:g5siyxDy:sym=",
             "=der:g5siyxDz:sym=", "=der:g5siyyDz:sym=", "=der:g5siyzDz:sym=",
             "=der:g5siz0D0:sym=", "=der:g5siz0Dx:sym=", "=der:g5siz0Dy:sym=",
             "=der:g5siz0Dz:sym=", "=der:g5sizxDx:sym=", "=der:g5sizxDy:sym=",
             "=der:g5sizxDz:sym=", "=der:g5sizyDy:sym=", "=der:g5sizyDz:sym="]
thrp_der2_gams = ["der2:g0DxDy", "der2:g0DxDz", "der2:g0DyDz"]
thrp_der3_gams = ["der3:g0D0DxDx", "der3:g0D0DyDy", "der3:g0D0DzDz",
                  "der3:g0DxDyDz", "der3:gxDxDyDy", "der3:gxDxDzDz",
                  "der3:gyDyDzDz"]

def momenta(max_qsq):
    moms = []
    for pz in range(-max_qsq, max_qsq+1):
        for py in range(-max_qsq, max_qsq+1):
            for px in range(-max_qsq, max_qsq+1):
                qsq = px**2+py**2+pz**2
                if qsq <= max_qsq:
                    moms.append([qsq, px, py, pz])
    for i in range(len(moms)):
        for j in range(i, len(moms)):
            if moms[j][0] <= moms[i][0]:
                swap = moms[i][:]
                moms[i][:] = moms[j][:]
                moms[j][:] = swap[:]    
    return moms

def get_data_type(fname):
    with h5py.File(fname, "r") as fp:
        # First of all there should be a lowest-level "arr" dataset
        names = []
        fp.visititems(lambda x,t: names.append(x) if 'arr' in x and type(t) is h5py.Dataset else None)
        assert names != [], " No 'arr' datasets found"
        if '=1=' in " ".join(names):
            return "meson"
        if 'fwd' in " ".join(names):
            return "proj-nucl"
        if '2-2' in " ".join(names):
            return "nucl"
        if 'twop' in " ".join(names):
            return "twop"
        if '=loc:1=' in " ".join(names):
            return "thrp"
        if 'der2:g0DxDy' in " ".join(names):
            return "thrp_der2"
        if 'der3:g0D0DxDx' in " ".join(names):
            return "thrp_der3"

def totxt_nproj(fname):
    with h5py.File(fname, "r") as fp:
        names = []
        fp.visititems(lambda x,t: names.append(x) if 'arr' in x and type(t) is h5py.Dataset else None)
        top = "/".join(names[0].split("/")[:2])
        msqs = list(sorted(set([n.split("/")[-2] for n in names])))
        max_qsq = max(msqs)
        max_qsq = int(re.match('msq0*([1-9][0-9]*)', max_qsq).groups()[0])
        data = {}
        for msq in msqs:
            for flv in nucl_flvs:
                for d in ["fwd", "bwd"]:                
                    grp = "%s/%s/%s/%s/%s" % (top, flv, '1-1', d, msq)
                    arr = np.array(fp[grp]['arr'])
                    mom = np.array(fp[grp]['mvec'])
                    data[msq,flv,d] = {"arr": arr, "mvec": mom}
        T = list(data.values())[0]['arr'].shape[0]
        momvecs = momenta(max_qsq)
        lines = []
        for t in range(T):
            for mv in momvecs:
                msq = "msq%04d" % mv[0]
                mom = data[msq,nucl_flvs[0],"fwd"]["mvec"]
                idx = mom.tolist().index(mv[1:])                                
                line = "  %2d %+d %+d %+d" % (t, mv[1], mv[2], mv[3])
                for flv in nucl_flvs:
                    for d in ["fwd", "bwd"]:
                        arr = data[msq,flv,d]["arr"]
                        mom = data[msq,flv,d]["mvec"]
                        assert idx == mom.tolist().index(mv[1:]), " momenta missmatch"
                        reim = arr[t,idx]
                        line += "   %+e %+e" % (reim.real, reim.imag)
                lines.append(line)
        buf = "\n".join(lines) + "\n"
        return buf

def totxt_nucl(fname):
    with h5py.File(fname, "r") as fp:
        names = []
        fp.visititems(lambda x,t: names.append(x) if 'arr' in x and type(t) is h5py.Dataset else None)
        top = "/".join(names[0].split("/")[:-4])
        msqs = list(sorted(set([n.split("/")[-2] for n in names])))
        max_qsq = max(msqs)
        max_qsq = int(re.match('msq0*([1-9][0-9]*)', max_qsq).groups()[0])
        data = {}
        for msq in msqs:
            for flv in nucl_flvs:
                for intrp in nucl_intrps:
                    grp = "%s/%s/%s/%s" % (top, flv, intrp, msq)
                    arr = np.array(fp[grp]['arr'])
                    mom = np.array(fp[grp]['mvec'])
                    data[msq,flv,intrp] = {"arr": arr, "mvec": mom}
        T = list(data.values())[0]['arr'].shape[0]
        momvecs = momenta(max_qsq)
        lines = []
        for t in range(T):
            for mv in momvecs:
                msq = "msq%04d" % mv[0]
                for flv in nucl_flvs:
                    for intrp in nucl_intrps:
                        arr = data[msq,flv,intrp]["arr"]
                        mom = data[msq,flv,intrp]["mvec"]
                        idx = mom.tolist().index(mv[1:])
                        reim = arr[t,idx,:,:]
                        for i in range(NS):                            
                            line = "  %2d %+d %+d %+d " % (t, mv[1], mv[2], mv[3])
                            for j in range(NS):
                                line += " %+e %+e " % (reim[i,j].real, reim[i,j].imag)
                            line += " %s %s" % (flv, intrp)
                            lines.append(line)
        buf = "\n".join(lines) + "\n"
        return buf

def totxt_meson(fname):
    with h5py.File(fname, "r") as fp:
        names = []
        fp.visititems(lambda x,t: names.append(x) if 'arr' in x and type(t) is h5py.Dataset else None)
        top = "/".join(names[0].split("/")[:-4])
        msqs = list(sorted(set([n.split("/")[-2] for n in names])))
        max_qsq = max(msqs)
        max_qsq = int(re.match('msq0*([1-9][0-9]*)', max_qsq).groups()[0])
        data = {}
        for msq in msqs:
            for gam in meson_gams:
                for flv in meson_flvs:
                    grp = "%s/%s/%s/%s" % (top, flv, gam, msq)
                    arr = np.array(fp[grp]['arr'])
                    mom = np.array(fp[grp]['mvec'])
                    data[msq,gam,flv] = {"arr": arr, "mvec": mom}
        T = list(data.values())[0]['arr'].shape[0]
        momvecs = momenta(max_qsq)
        lines = []
        for t in range(T):
            for mv in momvecs:
                msq = "msq%04d" % mv[0]
                for gam in meson_gams:
                    for flv in meson_flvs:
                        arr = data[msq,gam,flv]["arr"]
                        mom = data[msq,gam,flv]["mvec"]
                        idx = mom.tolist().index(mv[1:])
                        reim = arr[t,idx]
                        line = "  %2d %+d %+d %+d %+e %+e \t%s\t%s" % (t, mv[1], mv[2], mv[3],
                                                                       reim.real, reim.imag,
                                                                       gam, flv)
                        lines.append(line)
        buf = "\n".join(lines) + "\n"
        return buf

def totxt_twop(fname):
    with h5py.File(fname, "r") as fp:
        names = []
        fp.visititems(lambda x,t: names.append(x) if 'arr' in x and type(t) is h5py.Dataset else None)
        top = "/".join(names[0].split("/")[:-2])
        msqs = list(sorted(set([n.split("/")[-2] for n in names])))
        max_qsq = max(msqs)
        max_qsq = int(re.match('msq0*([1-9][0-9]*)', max_qsq).groups()[0])
        data = {}
        for msq in msqs:
            grp = "%s/%s" % (top, msq)
            arr = np.array(fp[grp]['arr'])
            mom = np.array(fp[grp]['mvec'])
            data[msq] = {"arr": arr, "mvec": mom}
        T = list(data.values())[0]['arr'].shape[0]
        momvecs = momenta(max_qsq)
        lines = []
        for t in range(T):
            for mv in momvecs:
                msq = "msq%04d" % mv[0]
                arr = data[msq]["arr"]
                mom = data[msq]["mvec"]
                idx = mom.tolist().index(mv[1:])
                reim = arr[t,idx]
                line = "  %2d %+d %+d %+d %+e %+e" % (t, mv[1], mv[2], mv[3],
                                                      reim.real, reim.imag)
                lines.append(line)
        buf = "\n".join(lines) + "\n"
        return buf

def totxt_thrp(fname):
    with h5py.File(fname, "r") as fp:
        names = []
        fp.visititems(lambda x,t: names.append(x) if 'arr' in x and type(t) is h5py.Dataset else None)
        top = "/".join(names[0].split("/")[:-3])
        msqs = list(sorted(set([n.split("/")[-2] for n in names])))
        max_qsq = max(msqs)
        max_qsq = int(re.match('msq0*([1-9][0-9]*)', max_qsq).groups()[0])
        data = {}
        for msq in msqs:
            for gam in thrp_gams:
                grp = "%s/%s/%s/" % (top, gam, msq)
                arr = np.array(fp[grp]['arr'])
                mom = np.array(fp[grp]['mvec'])
                data[msq,gam] = {"arr": arr, "mvec": mom}
        T = list(data.values())[0]['arr'].shape[0]
        momvecs = momenta(max_qsq)
        lines = []
        for t in range(T):
            for mv in momvecs:
                msq = "msq%04d" % mv[0]
                for gam in thrp_gams:
                    arr = data[msq,gam]["arr"]
                    mom = data[msq,gam]["mvec"]
                    idx = mom.tolist().index(mv[1:])
                    reim = arr[t,idx]
                    line = "  %2d %+d %+d %+d  %+e %+e %s" % (t, mv[1], mv[2], mv[3],
                                                              reim.real, reim.imag,
                                                              gam)
                    lines.append(line)
        buf = "\n".join(lines) + "\n"
        return buf

def totxt_thrp_der2(fname):
    with h5py.File(fname, "r") as fp:
        names = []
        fp.visititems(lambda x,t: names.append(x) if 'arr' in x and type(t) is h5py.Dataset else None)
        top = "/".join(names[0].split("/")[:-3])
        msqs = list(sorted(set([n.split("/")[-2] for n in names])))
        max_qsq = max(msqs)
        max_qsq = int(re.match('msq0*([1-9][0-9]*)', max_qsq).groups()[0])
        data = {}
        for msq in msqs:
            for gam in thrp_der2_gams:
                grp = "%s/%s/%s/" % (top, gam, msq)
                arr = np.array(fp[grp]['arr'])
                mom = np.array(fp[grp]['mvec'])
                data[msq,gam] = {"arr": arr, "mvec": mom}
        T = list(data.values())[0]['arr'].shape[0]
        momvecs = momenta(max_qsq)
        lines = []
        for t in range(T):
            for mv in momvecs:
                msq = "msq%04d" % mv[0]
                for gam in thrp_der2_gams:
                    arr = data[msq,gam]["arr"]
                    mom = data[msq,gam]["mvec"]
                    idx = mom.tolist().index(mv[1:])
                    reim = arr[t,idx]
                    line = "  %2d %+d %+d %+d  %+e %+e %s" % (t, mv[1], mv[2], mv[3],
                                                              reim.real, reim.imag,
                                                              gam)
                    lines.append(line)
        buf = "\n".join(lines) + "\n"
        return buf
    
def totxt_thrp_der3(fname):
    with h5py.File(fname, "r") as fp:
        names = []
        fp.visititems(lambda x,t: names.append(x) if 'arr' in x and type(t) is h5py.Dataset else None)
        top = "/".join(names[0].split("/")[:-3])
        msqs = list(sorted(set([n.split("/")[-2] for n in names])))
        max_qsq = max(msqs)
        max_qsq = int(re.match('msq0*([1-9][0-9]*)', max_qsq).groups()[0])
        data = {}
        for msq in msqs:
            for gam in thrp_der3_gams:
                grp = "%s/%s/%s/" % (top, gam, msq)
                arr = np.array(fp[grp]['arr'])
                mom = np.array(fp[grp]['mvec'])
                data[msq,gam] = {"arr": arr, "mvec": mom}
        T = list(data.values())[0]['arr'].shape[0]
        momvecs = momenta(max_qsq)
        lines = []
        for t in range(T):
            for mv in momvecs:
                msq = "msq%04d" % mv[0]
                for gam in thrp_der3_gams:
                    arr = data[msq,gam]["arr"]
                    mom = data[msq,gam]["mvec"]
                    idx = mom.tolist().index(mv[1:])
                    reim = arr[t,idx]
                    line = "  %2d %+d %+d %+d  %+e %+e %s" % (t, mv[1], mv[2], mv[3],
                                                              reim.real, reim.imag,
                                                              gam)
                    lines.append(line)
        buf = "\n".join(lines) + "\n"
        return buf
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("FNAME", type=str)
    parser.add_argument("-o", "--output", metavar="F", type=str, help="output file name (default: FNAME.txt)")
    args = parser.parse_args()    
    fname = args.FNAME
    output = args.output
    if output is None:
        output = fname + ".txt"

    ### What kind of file are we parsing?
    dt = get_data_type(fname)
    print(" %s file detected" % dt)
    if dt == "proj-nucl":
        buf = totxt_nproj(fname)
    if dt == "nucl":
        buf = totxt_nucl(fname)
    if dt == "meson":
        buf = totxt_meson(fname)
    if dt == "twop":
        buf = totxt_twop(fname)
    if dt == "thrp":
        buf = totxt_thrp(fname)
    if dt == "thrp_der2":
        buf = totxt_thrp_der2(fname)
    if dt == "thrp_der3":
        buf = totxt_thrp_der3(fname)

    with open(output, "w") as fp:
        print(" writing to %s" % output)
        fp.write(buf)
    return 0

if __name__ == "__main__":
    sys.exit(main())
