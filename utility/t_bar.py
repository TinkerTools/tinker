#!/usr/bin/env python
'''
Analysis of Tinker BAR data using pymbar.
Enter help for more info.
'''
# Requires pymbar: conda install -c conda-forge pymbar

HELP_MSG = '''
Quick start

  The output is printed as a table.
  [equ, dF] and [equ, sd_dF] are the estimated free energy difference and its uncertainty using equilibrated samples;
  [blockave_uncorr, sd_dF] is the standard error estimated using blocks in one file or using different files.
  Take caution when the difference between [equ, dF] and [uncorr, dF] is much larger than [equ, sd_dF(kcal/mol)].
'''

HELP_SUMMARY = '''
  SUMMARY
  --------
  dF_equ           free energy difference using equilibrated data (kcal/mol)
  dF_all           free energy difference using all data (kcal/mol)
  se_dF_block      standard error using block average (kcal/mol)
  se_dF_bar        standard error using BAR (kcal/mol)
  dF_fwd           Forward FEP (kcal/mol)
  dF_bwd           Backward FEP (kcal/mol)
  ddF_fwd_bwd      difference between fwd and bwd FEP (kcal/mol)
  dH               enthalpy difference using equilibrated data (kcal/mol)
  TdS              entropy difference using equilibrated data (kcal/mol)
  se_dH            standard error of dH (kcal/mol)
  sd_dE            sd of deltaE (kcal/mol)
  overlap_all      connectivity using all data
  overlap_equ      connectivity using equilibrated data
  nA               effective sample size
  nB

'''
HELP_MSG = '''
Details
  %s
  Full table
  --------
  start_A:         first frame of trajectory A  
  end_A:           last frame of trajectory A  
  g_A:             statistical inefficiency, correlation. 1 means the original samples are independent, or the correlation is not calculated.
  dF(kcal/mol):    BAR free energy difference in kcal/mol
  sd_dF(kcal/mol): standard error
  dF_fwd:          Forward FEP
  dF_bwd:          Backward FEP
  dF_fwd:          mean(deltaE) in forward FEP
  sd_dF_fwd:       sd(deltaE) in forward FEP
  p_A(traj_B):     probabiliy that samples in state A are from trajetory B
  overlap:         connectivity between two states; 1 means perfect overlap
  
  all:             using all samples. This assumes independent samples, so the free energy results will have larger variance and the uncertainty will be underestimated. 
  uncorr:          using uncorrelated samples. 
  equ:             using uncorrelated samples after equilibration. 

  block?xxx:       results for individual blocks
  blockave_xxx:    average and standard error calculated from blocks
  
  https://pymbar.readthedocs.io/en/master/mbar.html

'''%HELP_SUMMARY
# Email: francis {at} qubit-pharmaceuticals {dot} com

import numpy as np
import pandas as pd
import pymbar
import sys
import os
import argparse

def convert_1d(arr):
  '''
  convert nx2 array to nx1 array of the difference
  '''
  if len(arr.shape) == 1:
    arr1 = arr
  elif len(arr.shape) ==2 and arr.shape[1] == 2:
    arr1 = arr[:, 0] - arr[:, 1]
  else:
    arr1 = arr[:, 0]
  return arr1
  
def calc_eff_size(arr, equil=False):
  '''
  return t0, g
  '''
  arr1 = convert_1d(arr)
  if len(arr1) == 0:
    return 0, 1
  NBLOCK = 71
  nskip = max(1, len(arr1)//NBLOCK)
  t0 = 0
  indices = np.arange(len(arr1), dtype=int)
  if equil:
    [t0, g, Neff_max] = pymbar.timeseries.detectEquilibration(arr1, nskip=nskip)
    #indices = pymbar.timeseries.subsampleCorrelatedData(arr1[t0:], g=g)
  else:
    indices = pymbar.timeseries.subsampleCorrelatedData(arr1)
    g = len(arr1)/max(1, len(indices))
  return t0, g

def resample(arr, n, uniform=False):
  if uniform:
    idx = np.array((np.arange(n, dtype=int)*(len(arr)-1))//(n-1), dtype=int)
    return arr[idx]
  else:
    return arr[np.random.randint(0, len(arr), n, dtype=int)]
    
def subsample1(arr1, equil=False, corr=True, skip=0, last=None, n_aug=1, uniform=True):
  '''
  equil: discard equilibration samples
  corr: compute correlation
  skip: first sample
  last: last sample
  n_aug: multiply sample size by this factor
  '''
  kwargs = {'equil':equil}
  if skip is None:
    skip = 0
  if last is None:
    last = len(arr1)
  elif last < 0:
    last = last % len(arr1)
  if (skip-1)*skip<=0 and (last-1)*last<=0:
    skip = int(skip*len(arr1))
    last = int(last*len(arr1))
  last = min(last, len(arr1))

  arr1b = arr1[skip:last]
  t1, g1 = calc_eff_size(arr1b, **kwargs)
  if not corr:
    g1 = 1
  n1 = len(arr1b)*n_aug//max(1, g1)
  arr1c = resample(arr1b, n1, uniform=True)
  return arr1c, t1 + skip, last, g1

def subsample_arrays(arr1s, **kwargs):
  arr2s = []
  # t0, t1, g, n
  tgn = np.zeros((len(arr1s), 4))
  for i, arr1 in enumerate(arr1s):
    res = subsample1(arr1, **kwargs)
    arr2s.append(res[0])
    tgn[i, 0] = res[1]
    tgn[i, 1] = res[2]
    tgn[i, 2] = res[3]
    tgn[i, 3] = len(res[0])

  arr2 = np.concatenate(arr2s, axis=0)
  t1 = np.sum(tgn[:, 0])
  t2 = np.sum(tgn[:, 1])
  g1 = np.sum(tgn[:, 3] * tgn[:, 2])/max(1, np.sum(tgn[:, 3]))
  return arr2, t1, t2, g1

def subsample2(arr1, arr2, equil=False, corr=True, skip1=0, skip2=0, nmin=100):
  '''
  resample arrays based on their effective sample size
  the output arrays are proportional (but not equal) to the number of effetive samples

  return:
    arr1b: resampled arr1
    arr2b
    t1: starting frame number for arr1
    t2
    g1: effective sample size for arr1
    g2
    teff: effective sample size for resampled arrays
  '''
  kwargs = {'equil':equil}
  g1 = 1
  g2 = 1
  if equil:
    t1, g1 = calc_eff_size(arr1, **kwargs)
    t2, g2 = calc_eff_size(arr2, **kwargs)
    if not corr:
      g1 = 1
      g2 = 1
  else:
    t1 = skip1
    t2 = skip2
    if corr:
      _t1, g1 = calc_eff_size(arr1[skip1:], **kwargs)
      _t2, g2 = calc_eff_size(arr2[skip2:], **kwargs)
  teff = min(g1, g2)
  n1 = len(arr1) - t1
  n2 = len(arr2) - t2
  if corr:
    if g1 < g2:
      n1 = int(g2/g1*n1)
    elif g2 < g1:
      n2 = int(g1/g2*n2)
  _min = min(n1, n2)
  if _min < nmin:
    n1 = n1 * nmin//_min
    n2 = n2 * nmin//_min
    #print("Augmenting data %d > %d"%(_min, nmin))
  arr1b = resample(arr1[t1:], n1)
  arr2b = resample(arr2[t2:], n2)
  teff = n1/(len(arr1)-t1)*g1
  return arr1b, arr2b, t1, g1, t2, g2, teff

def subsample(arr, equil=False, corr=True):
  arr1 = convert_1d(arr)
  NBLOCK = 71
  nskip = max(1, len(arr1)//NBLOCK)
  t0 = 0
  indices = np.arange(len(arr1), dtype=int)
  if equil and corr:
    [t0, g, Neff_max] = pymbar.timeseries.detectEquilibration(arr1, nskip=nskip)
    indices = pymbar.timeseries.subsampleCorrelatedData(arr1[t0:], g=g)
  elif corr:
    indices = pymbar.timeseries.subsampleCorrelatedData(arr1)
  return [_+t0 for _ in indices]

def get_index_summary(indices):
  return '(%d,%d,%d)'%(indices[0], indices[-1], indices[1]-indices[0])

def get_index_summary2(idx1, idx2):
  msg = 'A:%s B:%s'%(get_index_summary(idx1), get_index_summary(idx2))
  return msg

def read_tinker_bars(inplist, temp=298, press=1.0, ibegin=0, iend=-1):
  u1s = []
  u2s = []
  temps = []
  for inp in inplist:
    res = read_tinker_bar(inp)
    if res is not None:
      res0 = res[0][ibegin:iend]
      res1 = res[1][ibegin:iend]
      #if len(res[0][ibegin:iend]) + len(res[1][ibegin:iend]) > 0:
      if len(res0) + len(res1) > 0:
        u1s.append(res0)
        u2s.append(res1)
        temps.append(res[2])
  if len(u1s) == 0:
    return None
  return u1s, u2s, temps

def read_tinker_bar(inp, temp=298, press=1.0):
  #beta = 1.0/(8.314*temp/4184)
  # bar * ang^3 in kcal/mol
  PV_UNIT_CONV = 1e5*1e-30/4184*6.02e23
  with open(inp, 'r') as fh:
    lines = fh.readlines()
    if len(lines) == 0:
      return
    w = lines[0].split()
    if len(w) == 0 or (not w[0].isdigit()):
      return
    n1 = int(w[0])
    temp1 = float(w[1]) if len(w) >= 2 else temp

    w = lines[n1+1].split()
    if len(w) == 0 or (not w[0].isdigit()):
      return
    n2 = int(w[0])
    temp2 = float(w[1]) if len(w) >= 2 else temp

    beta1 = 1.0/(8.314*temp1/4184)
    beta2 = 1.0/(8.314*temp2/4184)

    if len(lines) != n1+n2+2:
      print("nlines (%d) != %d + %d"%(len(lines), n1, n2))
      return

    arr1 = np.fromstring(''.join(lines[1:n1+1]), sep=' ').reshape(-1, 4)
    arr2 = np.fromstring(''.join(lines[n1+2:]), sep=' ').reshape(-1, 4)

    u1 = beta1*(arr1[:, 1:3] + press*arr1[:, 3:4]*PV_UNIT_CONV)
    u2 = beta2*(arr2[:, 1:3] + press*arr2[:, 3:4]*PV_UNIT_CONV)
    #return (np.concatenate((u1[idx1], u2[idx2]), axis=0).transpose(), [len(idx1), len(idx2)], msg)
    return u1, u2, temp1

def concat_arr(arr1, arr2, idx1, idx2):
  msg = get_index_summary2(idx1, idx2)
  return (np.concatenate((arr1[idx1], arr2[idx2]), axis=0).transpose(), [len(idx1), len(idx2)], msg)

def tinker_to_mbar(arr1, arr2, equil=False, corr=True):
  assert len(arr1.shape) == 2 and arr1.shape[1] == 2
  idx1 = subsample(arr1[:, 1] - arr1[:, 0], equil=equil, corr=corr)
  idx2 = subsample(arr2[:, 1] - arr2[:, 0], equil=equil, corr=corr)
  return concat_arr(arr1, arr2, idx1, idx2)

def arr_to_mbar(arr1, arr2):
  return (np.concatenate((arr1, arr2), axis=0).transpose(), [len(arr1), len(arr2)])

def exp_ave(ener, return_sd=True, ener_unit=1):
  if len(ener) == 0:
    return np.nan, 0
  emin = np.min(ener)
  eave = -np.log(np.mean(np.exp(-(ener-emin)))) + emin
  if return_sd:
    _sd = np.std(ener)
    return eave*ener_unit, _sd*ener_unit
  return eave*ener_unit

def calc_mbar(u_kn, N_k, teff=1, temp=298):
  #data = concat_arr(arr1, arr2, idx1, idx2)
  #cols = 'dF(kcal/mol) se_dF(kcal/mol) dH se_dH TdS se_TdS dF_fwd dF_bwd dE_fwd dE_bwd sd_dE_fwd sd_dE_bwd p_A(traj_B) p_B(traj_A) eig_overlap'.split()
  cols = 'dF(kcal/mol) se_dF(kcal/mol) dH TdS se_dH dF_fwd dF_bwd dE_fwd dE_bwd sd_dE_fwd sd_dE_bwd overlap'.split()
  if u_kn is None:
    return cols
  kt = temp*8.314/4184

  mbar = pymbar.MBAR(u_kn, N_k)
  #results = mbar.getFreeEnergyDifferences()
  results = mbar.computeEntropyAndEnthalpy()
  overlap = mbar.computeOverlap()
  #msg1 = get_index_summary(idx1)
  #msg2 = get_index_summary(idx2)
  es_fwd = (u_kn[1, 0:N_k[0]] - u_kn[0, 0:N_k[0]])
  es_bwd = (u_kn[0, N_k[0]:sum(N_k[0:2])] - u_kn[1, N_k[0]:sum(N_k[0:2])])
  de_fwd = np.mean(es_fwd)*kt
  de_bwd = np.mean(es_bwd)*kt

  fwd, sd_fwd = exp_ave(es_fwd, ener_unit=kt)
  bwd, sd_bwd = exp_ave(es_bwd, ener_unit=kt)
  ghs = [] # free energy, enthalpy, entropy
  for i in range(3):
    ghs.append(kt*results[i*2][0, 1])
    ghs.append(kt*results[i*2+1][0, 1]*np.sqrt(teff))
  #return ghs[0], ghs[1], ghs[2], ghs[3], ghs[4], ghs[5], fwd, -bwd, de_fwd, -de_bwd, sd_fwd, sd_bwd, overlap['matrix'][0, 1], overlap['matrix'][1, 0], overlap['scalar']
  return ghs[0], ghs[1], ghs[2], ghs[4], ghs[3], fwd, -bwd, de_fwd, -de_bwd, sd_fwd, sd_bwd, overlap['scalar']


def get_bar_res(data, ishift=0, **kwargs):
  '''
  data: [array1, array2, temp]
  return 1 row of data frame
  '''
  col0 = 'start_A end_A start_B end_B g_A g_B'.split()
  if data is None:
    return col0 + calc_mbar(None, None)
  opt1 = kwargs
  arr1, t1, t1b, g1 = subsample_arrays(data[0], **opt1)
  arr2, t2, t2b, g2 = subsample_arrays(data[1], **opt1)
  if len(arr1) + len(arr2) == 0:
    print("WARNING: empty array")
    return np.nan
  u_kn, N_k = arr_to_mbar(arr1, arr2)
  if len(u_kn) == 0:
    print("WARNING: empty array")
    return np.nan
  res = calc_mbar(u_kn, N_k, temp=data[2][0])
  return [_+ishift for _ in [t1, t1b, t2, t2b]] + [g1, g2] + list(res)


def compute_total_sd(df1):
  ''' compute sum on df columns
  sqrt(sum of squares) for 'sd_XXX', 'se_XXX'
  mininum for 'overlap*'
  '''
  arr1 = df1.sum(axis=0)
  for i, col in enumerate(df1.columns):
    if col.startswith('se_') or col.startswith('sd_'):
      arr1.iloc[i] = np.sqrt(np.sum(df1[col]**2.0))
    elif col.startswith('overlap') or col.endswith('overlap'):
      arr1.iloc[i] = df1[col].min()
    elif col.startswith('g_'):
      arr1.iloc[i] = df1[col].mean()
  return arr1

def compute_block_ave(df1):
  ''' compute average on df columns
  if df contains se_XXX and XXX, se_XXX will be calculated by the stderr of XXX
  '''
  arr1 = df1.mean(axis=0)
  nrow = len(df1.index)
  for i, col in enumerate(df1.columns):
    if col.startswith('se_') and col[3:] in df1.columns:
      arr1.iloc[i] = (df1[col[3:]].std())/np.sqrt(nrow)
  return arr1

def sum_df_list(dflist):
  ''' compute sum for list of dfs
  sqrt(sum of squares) for 'sd_XXX', 'se_XXX'
  mininum for 'overlap*' '*overlap'
  '''
  if len(dflist) == 0:
    return None
  elif len(dflist) == 1:
    return dflist[0]

  idx0 = dflist[0].index
  for df1 in dflist[1:]:
    idx0 = idx0 & df1.index
  df2list = []
  df3 = pd.DataFrame(columns=dflist[0].columns)
  for row in idx0:
    df2 = pd.concat([_.loc[[row], :] for _ in dflist], axis=0, ignore_index=True)
    df3.loc[row, :] = (compute_total_sd(df2))
  return df3

def check_data_size(data, ibegin, iend):
  data1 = [[], [], []]
  for i in range(len(data[0])):
    if len(data[0][i][ibegin:iend]) + len(data[1][i][ibegin:iend]) > 0:
      for j in range(len(data)):
        data1[j].append(data[j][i])
  return data1
def get_summary(df1, verbose=True):
  name_val = []
  name_val.append(('dF_equ', df1.loc['equ', 'dF(kcal/mol)']))
  name_val.append(('dF_all', df1.loc['all', 'dF(kcal/mol)']))
  if 'blockave_all' in df1.index:
    name_val.append(('se_dF_block', df1.loc['blockave_all', 'se_dF(kcal/mol)']))
  else:
    name_val.append(('se_dF_block', np.nan))
  nA = (df1.loc['equ', 'end_A']-df1.loc['equ', 'start_A'])//df1.loc['equ', 'g_A']
  nB = (df1.loc['equ', 'end_B']-df1.loc['equ', 'start_B'])//df1.loc['equ', 'g_B']

  sd_de_fwd = df1.loc['all', 'sd_dE_fwd']
  sd_de_bwd = df1.loc['all', 'sd_dE_bwd']

  name_val.append(('se_dF_bar', df1.loc['equ', 'se_dF(kcal/mol)']))
  name_val.append(('dF_fwd', df1.loc['all', 'dF_fwd']))
  name_val.append(('dF_bwd', df1.loc['all', 'dF_bwd']))
  name_val.append(('ddF_fwd_bwd', abs(df1.loc['all', 'dF_fwd']-df1.loc['all', 'dF_bwd'])))
  name_val.append(('dH', df1.loc['equ', 'dH']))
  name_val.append(('se_dH', df1.loc['equ', 'se_dH']))
  name_val.append(('sd_dE', np.sqrt(np.sum(np.array([sd_de_fwd, sd_de_bwd])**2.0))))
  name_val.append(('overlap_all', df1.loc['all', 'overlap']))
  name_val.append(('overlap_equ', df1.loc['equ', 'overlap']))
  name_val.append(('nA', nA))
  name_val.append(('nB', nB))
  names, vals = zip(*name_val)
  df2 = pd.DataFrame([vals], columns=names, index=['SUMMARY'])
  CHECK_THR = (('se_dF_block', 0.5, 1), ('se_dF_bar', 0.5, 1), ('ddF_fwd_bwd', 1.0, 1), ('sd_dE', 2.0, 1), ('overlap_all', 0.3, -1))
  name_desc = {}
  name_desc.update({'dF_equ':'Free energy difference (equilibrated) (kcal/mol)'})
  name_desc.update({'dF_all':'Free energy difference (all data) (kcal/mol) '})
  name_desc.update({'se_dF_block':'Standard error (block average) (kcal/mol) '})
  name_desc.update({'ddF_fwd_bwd':'Difference between forward and backward FEP (kcal/mol) '})
  name_desc.update({'sd_dE':'Standard deviation of energy difference (kcal/mol)'})
  name_desc.update({'overlap_all':'Connectivity between two states'})
  if verbose:
    for (_row, _thr, _sign) in CHECK_THR:
      if _row not in df2.columns:
        continue
      val = df2.loc['SUMMARY', _row] 
      if (val - _thr)*_sign > 0:
        _name = name_desc[_row] if _row in name_desc else _row
        if _sign > 0:
          _direct = 'above'
        else:
          _direct = 'below'
        print("WARNING: %s %s %s (%.3f vs. %.3f)"%(_name,  _direct, 'threshold', val, _thr))
  return df2

def calc_dg(flist, ibegin=0, iend=-1, nblocks=5, summary=True, verbose=True):
  NBLOCK = nblocks
  NMIN = 100
  names = 'all uncorr equ'.split()
  opts = [{'equil':False, 'corr':False}, 
          {'equil':False, 'corr':True},
          {'equil':True,  'corr':True},
          ]

  df_out = pd.DataFrame(columns=get_bar_res(None))
  namesb = names[:1]
  optsb = opts[:1]
  df_blocks = [pd.DataFrame(columns=get_bar_res(None)) for _ in namesb]
  n_name_blocks = 1
  data0 = read_tinker_bars(flist, ibegin=ibegin, iend=iend)
  if data0 is None:
    print("ERROR: empty data")
    return
  teff = 1
  opt1 = {'skip':0, 'last':None, 'uniform':True, "ishift":ibegin}
  for opt, name in zip(opts, names):
    opt1.update(opt)
    res = get_bar_res(data0, **opt1)
    if res is not None:
      df_out.loc[name, :] = res

  if len(data0[0]) == 1:
    # calculate results of blocks
    #n1 = len(data0[0][0][ibegin:iend])
    #n2 = len(data0[1][0][ibegin:iend])
    n1 = len(data0[0][0])
    n2 = len(data0[1][0])
    if min(n1, n2) >= NBLOCK*NMIN:
      #bn1 = n1//NBLOCK
      #bn2 = n2//NBLOCK
      iblk = -1
      for opt, name1 in zip(optsb, namesb):
        iblk += 1
        for i in range(NBLOCK):
          name = "block%d%s"%(i+1, name1)
          opt1.update({'skip':i/NBLOCK, 'last':(i+1)/NBLOCK})
          opt1.update(opt)
          res = get_bar_res(data0, **opt1)
          if res is not None:
            #df_out.loc[name, :] = res
            df_blocks[iblk].loc[name, :] = res
        break
  elif len(data0[0]) > 1:
    #opt1 = {'skip':ibegin, 'last':iend, 'uniform':True}
    opt1 = {'skip':0, 'last':None, 'uniform':True, "ishift":ibegin}
    iblk = -1
    for opt, name1 in zip(optsb, namesb):
      iblk += 1
    # calculate results of files
      for i in range(len(data0[0])):
        name = "file%d%s"%(i+1, name1)
        opt1.update(opt)
        data1 = [data0[0][i:i+1], data0[1][i:i+1], data0[2][i:i+1]]
        res = get_bar_res(data1, **opt1)
        if res is not None:
          df_blocks[iblk].loc[name, :] = res

  for name, df2 in zip(namesb, df_blocks):
    if len(df2.index) > 1:
      df_out.loc['blockave_%s'%name, :] = compute_block_ave(df2)
  df_out = pd.concat([df_out]+df_blocks)
  if verbose:
    print(df_out.to_string(max_rows=len(df_out.index), max_cols=len(df_out.columns)))
  df_sum = get_summary(df_out)
  if summary:
    print()
    print(df_sum.to_string(max_rows=len(df_sum.index), max_cols=len(df_sum.columns)))
  return df_out, df_sum, (namesb, df_blocks)

def add_line(inp):
  return '\n'+'*'*8+' '+inp+' '+'*'*8+'\n'
def list_to_string(alist, maxl=12):
  if len(alist) == 0:
    return '[]'
  def trim(s):
    if len(s) > maxl:
      return '..'+s[-maxl:]
    else:
      return s

  return '[%s ... %s]'%(trim(str(alist[0])), trim(str(alist[-1])))

def calc_dg_summary(flist, seperate=False, outfile=None, **kwargs):
  '''
  wrapper for calc_dg

  if seperate, analyze each file independently
  otherwise combine all files
  '''
  if seperate and len(flist) > 1:
    dflist = []
    df1list = []
    df_block_list = []
    names_block = []
    if outfile is None:
      fout = sys.__stdout__
    else:
      fout = open(outfile, 'w')
    for f1 in flist:
      print(add_line('Analysis of %s'%f1))
      df1, df1sum, data_block = calc_dg([f1], **kwargs)
      dflist.append(df1sum)
      df1list.append(df1)
      df_block_list.append(data_block[1])
      names_block = data_block[0]
    df_blocks = [sum_df_list(_) for _ in zip(*df_block_list)]
    print(add_line('Sum of %d files %s'%(len(flist), list_to_string(flist))))
    df3 = sum_df_list(df1list)
    for name, df_block in zip(names_block, df_blocks):
        df3.loc['blockave_%s'%name, :] = compute_block_ave(df_block)
    print(df3.to_string())
    df2 = pd.concat(dflist, axis=0)
    #df2 = pd.DataFrame(df2.to_numpy(), index=flist, columns=df2.columns)
    df2.set_axis(flist, axis=0, inplace=True)
    #df2.reindex(flist)
    #df2.loc['Total', :] = compute_total_sd(df2)
    df2.loc['Total', :] = get_summary(df3, verbose=False).iloc[0, :].to_numpy()
    for col in 'nA nB'.split():
      if col in (df2.columns & df3.columns):
        df2.loc['Total', col] = df2.iloc[:-1, col].sum()

    print(add_line('SUMMARY'))
    print(df2.to_string(), file=fout)
    if outfile is not None:
      fout.close()
  else:
    df1, df1sum, _ = calc_dg(flist, **kwargs)

def main():
  parser = argparse.ArgumentParser(prog='t_bar.py', usage='%(prog)s [options]', description=__doc__)
  parser.add_argument('barfiles', nargs='+', metavar='flist', help='list of Tinker .bar file(s)')
  parser.add_argument('-b', nargs='?', metavar='begin', type=int, default=0, required=False, help='first frame')
  parser.add_argument('-e', nargs='?', metavar='end', type=int, default=None, required=False, help='last frame')
  parser.add_argument('-nb', nargs='?', metavar='nblocks', type=int, default=5, required=False, help='number of blocks')
  parser.add_argument('-t', nargs='?', metavar='type', default='tinker', required=False, help='file type')
  parser.add_argument('-s',                           action='store_true', required=False, help='analyze each .bar file seperately')
  parser.add_argument('-o', nargs='?', metavar='output',                      required=False, help='output file if each .bar is analyzed seperately')
  parser.add_argument('-q',                             action='store_true', default=False, required=False, help='quiet mode')

  v = (parser.parse_args(sys.argv[1:]))
  if 'help' in v.barfiles:
    print(HELP_MSG)
    return
  if v.t != 'tinker':
    print("Only tinker format is supported")
    return
  #calc_dg(v.barfiles, ibegin=v.b, iend=v.e, nblocks=v.nb, summary=True, verbose=(not v.q))
  calc_dg_summary(v.barfiles, seperate=v.s, outfile=v.o, ibegin=v.b, iend=v.e, nblocks=v.nb, summary=True, verbose=(not v.q))
  
if __name__ == '__main__':
  main()
