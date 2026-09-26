#!/usr/bin/env python3
from __future__ import annotations
from dataclasses import dataclass
from itertools import combinations
from collections import defaultdict
import numpy as np
import scipy.linalg
import scipy.sparse as sp
from spin_fci import excitation_generator, annihilate, create

@dataclass
class Entry:
    pool_index:int; rank:int; group_key:tuple; source_indices:tuple
    coeffs:np.ndarray; specs:tuple; matrix:sp.csr_matrix
    leakage_ratio:float; commutator_ratio:float; fro_norm:float
    def as_dict(self):
        return {'pool_index':int(self.pool_index),'rank':int(self.rank),
        'group_key':[list(x) if isinstance(x,tuple) else x for x in self.group_key],
        'source_indices':[int(x) for x in self.source_indices],
        'coeffs':[float(x) for x in self.coeffs],
        'specs':[{'remove':[int(x) for x in r],'add':[int(x) for x in a]} for r,a in self.specs],
        'support_size':int(np.sum(np.abs(self.coeffs)>1e-12)),
        'leakage_ratio':float(self.leakage_ratio),'commutator_ratio':float(self.commutator_ratio),
        'fro_norm':float(self.fro_norm)}

def _canon(v,tol=1e-12):
    v=np.asarray(v,float).copy(); v[np.abs(v)<tol]=0.; n=np.linalg.norm(v)
    if n==0:return v
    v/=n; nz=np.flatnonzero(np.abs(v)>tol)
    if len(nz):
        p=nz[np.argmax(np.abs(v[nz]))]
        if v[p]<0:v*=-1
    targets=(0.,1.,-1.,1/np.sqrt(2.),-1/np.sqrt(2.))
    for i,x in enumerate(v):
        for t in targets:
            if abs(x-t)<5e-10:v[i]=t;break
    n=np.linalg.norm(v)
    return v/n if n else v

def _rank(cols,tol=1e-10):
    return 0 if not cols else int(np.linalg.matrix_rank(np.column_stack(cols),tol=tol))

def sparse_null_basis(L,rcond=1e-11,zero_tol=1e-11,max_subset=12):
    L=np.asarray(L); m=L.shape[1]; N=scipy.linalg.null_space(L,rcond=rcond); dim=N.shape[1]
    if dim==0:return []
    out=[]
    if m<=max_subset:
        for s in range(1,m+1):
            for sub in combinations(range(m),s):
                Ns=scipy.linalg.null_space(L[:,sub],rcond=rcond)
                for j in range(Ns.shape[1]):
                    c=np.zeros(m); c[list(sub)]=np.real_if_close(Ns[:,j]).astype(float); c=_canon(c,zero_tol)
                    if np.linalg.norm(c)==0:continue
                    if np.linalg.norm(L@c)>1e-8*max(1.,np.linalg.norm(L)):continue
                    if _rank(out+[c])>_rank(out):
                        out.append(c)
                        if len(out)==dim:return out
    for j in range(dim):
        c=_canon(np.real_if_close(N[:,j]).astype(float),zero_tol)
        if _rank(out+[c])>_rank(out):out.append(c)
        if len(out)==dim:return out
    raise RuntimeError('could not construct deterministic null basis')

def _orientation(rem,add,norb):
    R=tuple(sorted(x%norb for x in rem)); A=tuple(sorted(x%norb for x in add))
    if R<A:return tuple(rem),tuple(add),R,A
    if R>A:return tuple(add),tuple(rem),A,R
    return (tuple(rem),tuple(add),R,A) if tuple(rem)<tuple(add) else (tuple(add),tuple(rem),A,R)

def generalized_specs(norb,rank):
    if rank==1:
        specs=[]; groups=defaultdict(list)
        for p,q in combinations(range(norb),2):
            for spin in (0,1):
                rem=(spin*norb+p,); add=(spin*norb+q,); i=len(specs); specs.append((rem,add)); groups[(1,(p,),(q,))].append(i)
        return specs,groups
    if rank!=2:raise ValueError('rank must be 1 or 2')
    pairs=list(combinations(range(2*norb),2)); specs=[]; groups=defaultdict(list); seen=set()
    for r0 in pairs:
        nr=sum(x<norb for x in r0)
        for a0 in pairs:
            if set(r0)&set(a0) or sum(x<norb for x in a0)!=nr:continue
            rem,add,R,A=_orientation(r0,a0,norb); key=(rem,add)
            if key in seen:continue
            seen.add(key); i=len(specs); specs.append((rem,add)); groups[(2,R,A)].append(i)
    return specs,groups

def _local_splus(m):
    dim=1<<(2*m); rows=[];cols=[];vals=[]
    for d in range(dim):
        for p in range(m):
            x,s1=annihilate(d,m+p)
            if not s1:continue
            x,s2=create(x,p)
            if not s2:continue
            rows.append(x);cols.append(d);vals.append(s1*s2)
    return sp.csr_matrix((vals,(rows,cols)),shape=(dim,dim),dtype=float)

def _local_coeffs(specs,norb,rcond=1e-11):
    spat=sorted({x%norb for r,a in specs for x in r+a}); mp={o:i for i,o in enumerate(spat)}; m=len(spat)
    def lm(x):return (x//norb)*m+mp[x%norb]
    loc=[(tuple(lm(x) for x in r),tuple(lm(x) for x in a)) for r,a in specs]
    bas=np.arange(1<<(2*m),dtype=np.int64); Sp=_local_splus(m); Cs=[]
    for r,a in loc:
        A=excitation_generator(bas,r,a).tocsr(); Cs.append((Sp@A-A@Sp).tocsr())
    n=len(Cs); G=np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):G[i,j]=G[j,i]=float(Cs[i].multiply(Cs[j]).sum())
    return sparse_null_basis(G,rcond=rcond)

def build_pool(basis,U,splus,norb,nalpha,nbeta,max_rank=2):
    # 18q used the global doublet-leakage construction; 24q used the
    # equivalent faster local-Fock [S_+,B]=0 construction.
    use_fast=norb>9; ms=.5*(nalpha-nbeta)
    S2=(splus.conj().T@splus).tocsr()+ms*(ms+1)*sp.eye(len(basis),format='csr')
    V=None if use_fast else scipy.linalg.orth(splus.conj().T.toarray(),rcond=1e-12)
    out=[]
    for rank in range(1,max_rank+1):
        specs,groups=generalized_specs(norb,rank)
        cache={}
        def A(k):
            if k not in cache:cache[k]=excitation_generator(basis,*specs[k]).tocsr()
            return cache[k]
        for gkey,idxs in sorted(groups.items(),key=lambda kv:kv[0]):
            if use_fast:
                coeffs=_local_coeffs([specs[k] for k in idxs],norb)
            else:
                blocks=[V.conj().T@(A(k)@U) for k in idxs]; m=len(blocks); G=np.zeros((m,m))
                for i in range(m):
                    for j in range(i,m):G[i,j]=G[j,i]=float(np.vdot(blocks[i],blocks[j]).real)
                coeffs=sparse_null_basis(G)
            for c in coeffs:
                B=sp.csr_matrix((len(basis),len(basis)),dtype=float)
                for cc,k in zip(c,idxs):
                    if abs(cc)>1e-13:B=B+float(cc)*A(k)
                B.eliminate_zeros(); fro=float(sp.linalg.norm(B))
                if fro<=1e-12:continue
                comm=float(sp.linalg.norm(S2@B-B@S2))/max(fro,1e-30)
                leak=0. if use_fast else float(np.linalg.norm(V.conj().T@(B@U)))/max(float(np.linalg.norm(B@U)),1e-30)
                out.append(Entry(len(out),rank,gkey,tuple(idxs),np.asarray(c,float),tuple(specs[k] for k in idxs),B,leak,comm,fro))
    return out

def entry_from_dict(d,basis):
    specs=tuple((tuple(x['remove']),tuple(x['add'])) for x in d['specs']); coeff=np.asarray(d['coeffs'],float)
    B=sp.csr_matrix((len(basis),len(basis)),dtype=float)
    for c,s in zip(coeff,specs):
        if abs(c)>1e-13:B=B+float(c)*excitation_generator(basis,*s)
    return Entry(int(d['pool_index']),int(d['rank']),tuple(d.get('group_key',())),tuple(d.get('source_indices',())),coeff,specs,B.tocsr(),float(d.get('leakage_ratio',0)),float(d.get('commutator_ratio',0)),float(d.get('fro_norm',sp.linalg.norm(B))))

def fermionic_terms_for_entry(d):
    coeff=np.asarray(d['coeffs'],float); specs=d['specs']; terms={}
    for c,s in zip(coeff,specs):
        if abs(c)<=1e-13:continue
        rem=tuple(s['remove']); add=tuple(s['add'])
        fwd=' '.join([f'+_{p}' for p in add]+[f'-_{q}' for q in reversed(rem)])
        rev=' '.join([f'+_{q}' for q in rem]+[f'-_{p}' for p in reversed(add)])
        terms[fwd]=terms.get(fwd,0)+float(c); terms[rev]=terms.get(rev,0)-float(c)
    return {k:v for k,v in terms.items() if abs(v)>1e-13}
