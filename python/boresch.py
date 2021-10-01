#!/usr/bin/env python
'''
Get Tinker keywords for Boresch restraint

> get_rot_rest.py sample.xyz sample.key

Read any traj format and "ligand keyword" from Tinker key file
Print 6 distance/angle/torsion restraints and the standard state correction

'''

import numpy as np
import mdtraj as md
import sys
import os
import time
import scipy.optimize
from scipy.spatial import distance_matrix
from scipy.cluster.hierarchy import fcluster, leaders, linkage
from scipy.spatial.distance import pdist
from collections import Counter, defaultdict
from openbabel import openbabel
from openbabel import pybel

def get_adjlist(top):
    adjlist = defaultdict(list)
    for bond in top.bonds:
        a1, a2 = bond
        if not (a1.name.startswith('H') or a2.name.startswith('H') ):
            adjlist[a1.index].append(a2.index)
            adjlist[a2.index].append(a1.index)
    return adjlist

def get_tab_branch(top):
    '''
    Given mdtraj Trajectory.Topology, 
    returns a dictionary for the number of connected heavy atoms
    '''
    i = 0
    nbr = np.zeros(top.n_atoms, dtype=np.int)
    for bond in top.bonds:
        i += 1
        a1, a2 = bond
        if not (a1.name.startswith('H') or a2.name.startswith('H') ):
            nbr[a1.index] += 1
            nbr[a2.index] += 1
    return nbr

def error_exit(msg):
    print("ERROR:", msg)
    sys.exit(1)
def write_tinker_idx(idxs):
    idx_out = []
    rs = []
    for i0 in sorted(idxs):
        if len(rs) == 0 or rs[-1][1]+1 < i0:
            rs.append([i0, i0])
        else: 
            rs[-1][1] = i0
    for r in rs:
        if r[0] == r[1]:
            idx_out.append(r[0])
        elif r[0] == r[1] - 1:
            idx_out.append(r[0])
            idx_out.append(r[1])
        else:
            idx_out.append(-r[0])
            idx_out.append(r[1])
    return idx_out
        
def read_tinker_idx(args):
    ''' Read tinker indices

    Example:
    ['5'] -> {5}
    ['-5', '7', '10', '-12', '15'] -> {5, 6, 7, 10, 12, 13, 14, 15}

    args: list of strings of integers
    return: a set of indices

    '''
    _range = []
    idxs = []
    for a in args:
        n = int(a)
        if n < 0 and len(_range) == 0:
            _range.append(-n)
        if n > 0:
            if len(_range) == 1:
                idxs.extend(list(range(_range.pop(), n+1)))
            else:
                idxs.append(n)
    return set(idxs)

def read_ligidx(fkey):
    '''Read tinker key file and return indices of the ligand
    '''
    ligand_idx = set() # ligand indices
    with open(fkey, 'r') as fh:
        for line in fh:
            w = line.split()
            if line.lower().startswith('ligand') and len(w) >= 2:
                idx_str = line.replace(',', ' ').split()[1:]
                ligand_idx |= read_tinker_idx(idx_str)
    return ligand_idx


def target_disp(weights, coord, alpha=0.1):
    '''target function for center of mass displacement and number of atoms in the group
    '''
    assert len(weights) == coord.shape[0]
    wts = np.array(weights).reshape((-1, 1))
    wts = np.maximum(wts, 0)
    assert sum(wts) > 0
    wts *= 1.0/np.mean(wts)

    loss = np.sum(np.abs(np.sum(wts*coord, axis=0)))
    loss += alpha*np.sum(np.abs(np.abs(wts*coord)))

    for m in np.arange(2, 20):
        mask0 = (wts > m)
        loss += alpha*np.sum((np.abs((wts*mask0))))

    return loss

def calc_coord(traj, idxs):
    vals = np.zeros(traj.n_frames)
    idxs = np.array(idxs)
    NM_IN_ANG = 10.0
    RAD_IN_DEG = 180.0/np.pi
    if len(idxs) == 2:
        vals = md.compute_distances(traj, idxs.reshape(1, -1))
        vals *= NM_IN_ANG
    elif len(idxs) == 3:
        vals = md.compute_angles(traj, idxs.reshape(1, -1))
        vals *= RAD_IN_DEG
    elif len(idxs) == 4:
        vals = md.compute_dihedrals(traj, idxs.reshape(1, -1))
        vals *= RAD_IN_DEG
    return vals[:, 0]

def calc_idx_ortho(traj, i1, i2, idx3, nbr=None, method='long'):
    '''
    find the index that gives the largest ortho vector

    nbr: list of nr of branched atoms
    '''
    DELTA_R = 0.2
    DELTA_COS = 0.3

    if nbr is None:
        nbr = defaultdict(int)
    t = traj
    a1 = t.xyz[0, [i1], :] - t.xyz[0, [i2], :]
    u1 = a1/np.linalg.norm(a1, axis=1)

    vec2 = t.xyz[0, idx3, :] - t.xyz[0, [i2], :]
    d2 = np.linalg.norm(vec2, axis=1)
    d2p = np.abs(np.sum(vec2 * u1,axis=1))
    d2o = np.sqrt(d2**2.0 - d2p**2.0)
    d2cos = d2p / np.maximum(1e-5, d2)
    if method == 'short':
        # large angle (~90 deg), short distance
        d2on = np.array(list(zip([nbr[_]<1 for _ in idx3], -d2cos//DELTA_COS, -d2o)) , dtype=[('n', np.int),('c', np.float), ('r', np.float)])
    else:
        # large ortho vector
        d2on = np.array(list(zip([nbr[_] for _ in idx3], d2o)), dtype=[('n', np.int), ('r', np.float)])
    i3 = idx3[np.argsort(d2on)[-1]]
    return i3

def write_xyz_from_md(traj, idx):
    outp = '%d\nExtracted from MD\n'%(len(idx))
    for n, i in enumerate(idx):
        atom = traj.topology.atom(i)
        xyz = traj.xyz[0, i, :]*10
        outp += '%5s %12.6f %12.6f %12.6f\n'%(atom.element.symbol, xyz[0], xyz[1], xyz[2])
    return outp

def get_rotlist(traj, ligidx):
    ftype = 'xyz'
    ligidx_sort = sorted(ligidx)
    outp_xyz = write_xyz_from_md(traj, ligidx_sort)
    mymols = list([pybel.readstring(ftype, outp_xyz)])
    mymol = mymols[0]
    iter_bond = openbabel.OBMolBondIter(mymol.OBMol)
    rotlist = []
    for bond in iter_bond:
        if bond.IsRotor():
            i1, i2 = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            rotlist.append((ligidx_sort[i1-1], ligidx_sort[i2-1])) 
    return rotlist

def get_idx_pocket(traj, idx1, idx2, rcutoff=0.5):
    R_CLUSTER = rcutoff
    t = traj
    nbr = get_tab_branch(t.topology)
    sorted_i1, sorted_i2, dist1, dist2 = calc_idx_interface(traj, idx1, idx2, rcutoff=rcutoff)

    if len(sorted_i1) == 0:
        error_exit('No interface atom found within %.3f nm'%rcutoff)

    prot_idx = sorted_i1
    lig_idx = sorted_i2
    prot_dist = pdist(t.xyz[0, prot_idx, :]) 
    prot_distmat = distance_matrix(t.xyz[0, prot_idx, :],t.xyz[0, prot_idx, :]) 
    Z = linkage(prot_dist, method='average')

    clst_idx = fcluster(Z, R_CLUSTER, criterion='distance')
    clst_size = Counter(clst_idx)
    clst_size_list = sorted(clst_size.items(), key = lambda t:t[1])
    iclstm = clst_size_list[-1][0]

    flag_max = (clst_idx == iclstm)
    iprot = prot_idx[flag_max][0]

    lig_dist = np.linalg.norm(t.xyz[0, lig_idx, :] - t.xyz[0, [iprot], :], axis=1)
    dist_nr2 = np.array(list(zip(lig_dist, [nbr[_]<=1 for _ in lig_idx])), dtype=[('r', np.float), ('n', np.int)])
    ord_nr2 = np.argsort(dist_nr2, order=('n', 'r'))
    ilig = sorted_i2[ord_nr2[0]]

    iprot2 = calc_idx_ortho(traj, ilig, iprot, prot_idx, nbr)
    iprot3 = calc_idx_ortho(traj, iprot, iprot2, prot_idx[prot_idx != iprot], nbr)

    ilig2 = calc_idx_ortho(traj, iprot, ilig, lig_idx[lig_idx != ilig], nbr, method='short')
    ilig3 = calc_idx_ortho(traj, ilig, ilig2, lig_idx[(lig_idx != ilig)*(lig_idx != ilig2)], nbr, method='short')

    #                             -r-
    #                        -a-      -a-
    # iprot3 ... iprot2 -t- iprot -t- ilig -t- ilig2 ... ilig3
    int_idxs = [[iprot, ilig]]
    int_idxs.extend([[iprot2, iprot, ilig], [iprot, ilig, ilig2]])
    int_idxs.extend([[iprot3, iprot2, iprot, ilig], [iprot2, iprot, ilig, ilig2], [iprot, ilig, ilig2, ilig3]])
    #r0s = [calc_int(t.xyz[0, _, :])[0] for _ in int_idxs]
    #print(r0s)
    r0s = [calc_coord(t, _)[0] for _ in int_idxs]

    return int_idxs, r0s

def calc_idx_interface(traj, idx1, idx2, rcutoff=0.6):
    t = traj
    distmat = distance_matrix(t.xyz[0, idx1, :], t.xyz[0, idx2, :])
    dist2 = np.min(distmat, axis=0)
    dist1 = np.min(distmat, axis=1)

    ord1 = np.argsort(dist1)
    ord2 = np.argsort(dist2)

    n1 = np.sum(dist1 <= rcutoff)
    n2 = np.sum(dist2 <= rcutoff)

    return idx1[ord1[:n1]], idx2[ord2[:n2]], dist1[ord1[:n1]], dist2[ord2[:n2]]

def write_rest(int_idxs, r0s, k0s, fmt='tinker'):
    rest_name = ['', '', 'restrain-distance', 'restrain-angle', 'restrain-torsion']
    outp = ''
    for ridx, r0, k0 in zip(int_idxs, r0s, k0s):
        nat = len(ridx)
        if nat >= len(rest_name):
            continue
        if fmt == 'tinker':
            outp += '%s %s %.6f %.6f %.6f\n'%(rest_name[nat], ' '.join('%4d'%(_+1) for _ in ridx), k0, r0, r0)
        else:
            print("Format %s not supported"%fmt)
            return ''
    return outp

def get_rottors(traj, idx1):
    t = traj
    rotlist = get_rotlist(t, idx1)
    tors = []
    adjlist = get_adjlist(traj.topology)
    for bond in rotlist:
        a2 = min(bond)
        a3 = max(bond)
        a1s = [_ for _ in adjlist[a2] if _ != a3]
        a4s = [_ for _ in adjlist[a3] if _ != a2]
        if len(a1s) * len(a4s) == 0:
            error_exit("Cannot find heavy atom connected to rotable bond %d %d"%(a2, a3))
        tors.append([a1s[0], a2, a3, a4s[0]])
    r0s = [calc_coord(t, _)[0] for _ in tors]
    k0s = np.zeros_like(r0s) + 0.01
    print(write_rest(tors, r0s, k0s))
    return tors

def find_rotation_rest(fxyz, ligidx0, atomnames='CA', rcutoff=0.6, alpha=1.0, rest_rotbond=False, rotbond_only=False):
    try:
        t = md.load_arc(fxyz)
    except IOError:
        t = md.load(fxyz)
    if rest_rotbond:
        get_rottors(t, ligidx0)
    if rotbond_only:
        return

    nbr = get_tab_branch(t.topology)

    ligidx = np.array(sorted(list(set(ligidx0) - set(t.topology.select('name H')))))
    if len(ligidx) == 0:
        error_exit("No ligand heavy atoms found")

    protidx0 = t.topology.select('name %s'%(atomnames))
    protidx0 = np.array(sorted(set(protidx0) - set(ligidx)))

    if len(protidx0) == 0:
        error_exit("No ligand atoms within %.3f nm of protein %s atoms"%(rcutoff, atomnames))

    int_idxs, r0s = get_idx_pocket(t, protidx0, ligidx)
    RAD_IN_DEG = 180/np.pi
    k0s = np.zeros_like(r0s) + 10.0/(RAD_IN_DEG)**2.0
    k0s[0] = 10.0
    print(write_rest(int_idxs, r0s, k0s))
    RT = 8.314 * 298 / 4184
    # https://doi.org/10.1021/jp0217839
    # Eq. (14)
    assert len(k0s) == 6
    dgrest = RT*np.log(1662*8*np.pi**2.0*np.sqrt(np.prod(k0s)*RAD_IN_DEG**10.0) \
                        /(r0s[0]**2.0*np.sin(r0s[1]/RAD_IN_DEG)*np.sin(r0s[2]/RAD_IN_DEG)*(2*np.pi*RT)**3))

    print("#dGrest(kcal/mol) %.5f"%dgrest)



def find_grp_idx(fxyz, ligidx0, atomnames='CA', rcutoff=1.2, alpha=1.0):
    try:
        t = md.load_arc(fxyz)
    except IOError:
        t = md.load(fxyz)
    ligidx = np.array(sorted(list(set(ligidx0) - set(t.topology.select('name H')))))
    if len(ligidx) == 0:
        error_exit("No ligand heavy atoms found")

    protidx0 = t.topology.select('name %s'%(atomnames))
    protidx0 = np.array(sorted(set(protidx0) - set(ligidx)))

    distmat = distance_matrix(t.xyz[0, ligidx, :], t.xyz[0, protidx0, :])
    distpro = np.min(distmat, axis=0)
    protidx = protidx0[distpro <= rcutoff]
    if len(protidx) == 0:
        return

    imindist = np.argmin(distmat)
    iminlig = ligidx[imindist // len(protidx0)]

    com_lig = np.mean(t.xyz[0, list(ligidx), :], axis=0)
    com_lig = t.xyz[0, [iminlig], :]


    wts0 = np.ones(len(protidx))
    xyz_prot = t.xyz[0, protidx, :]
    xyz1_prot = xyz_prot - com_lig.reshape((1, -1))
    res = scipy.optimize.minimize(target_disp, wts0, args=(xyz1_prot, alpha))
    wts1 = np.array(res.x)
    wts1 = np.maximum(0, wts1)
    wts1 *= 1.0/np.sum(wts1)
    wtm = np.mean(wts1[wts1 > 0])
    mask1 = wts1 > 0.4*wtm
    #print(wts1)
    #print(np.mean(wts1))
    #print(mask1)
    com_p1 = (np.mean(t.xyz[0, protidx[mask1], :], axis=0)).reshape((1, -1))

    xyz_lig = t.xyz[0, ligidx, :]
    dists = np.linalg.norm(xyz_lig - com_p1, axis=1)
    #print('DIST', dists)
    imin = np.argmin(dists)

    dcom = np.linalg.norm(com_lig - com_p1)
    dmin = np.linalg.norm(xyz_lig[imin] - com_p1)


    idx_tinker = write_tinker_idx([_+1 for _ in protidx[mask1]])
    #print('#', (' '.join(['%5d'%(_) for _ in protidx[mask1]])))
    outp = ''
    sgrp = ''
    for n in idx_tinker:
        sgrp += ' %5d'%n
        if len(sgrp) > 50 and n > 0:
            #print('group 1 %s'%sgrp)
            outp += ('group 1 %s\n'%sgrp)
            sgrp = ''
            
    outp += ('group 2 %s\n'%(' '.join(['%5d'%(_) for _ in [ligidx[imin] + 1]])))
    outp += ("#r_0=%.3f"%(dmin*10))
    #print('group 2 %s'%(' '.join(['%5d'%(_) for _ in [ligidx[imin] + 1]])))
    #print("#r_0=%.3f"%(dmin*10))
    return dmin*10, outp


def main():
    if len(sys.argv) <= 2:
        print(__doc__)
        return 
    fxyz = sys.argv[1]
    fkey = sys.argv[2]
    ligidx = [_-1 for _ in sorted(read_ligidx(fkey))]
    args = sys.argv[3:]
    if len(ligidx) == 0:
        error_exit("ligand keyword not found")
    #rest_rotbond = ( len(sys.argv)>3 and sys.argv[3] == 'rotatable')
    rest_rotbond = 'rotatable' in args
    no_orient = 'only' in args
    #find_rotation_rest(fxyz, ligidx, "CA C N")
    find_rotation_rest(fxyz, ligidx, "CA C N", rest_rotbond=rest_rotbond, rotbond_only=no_orient)

main()

