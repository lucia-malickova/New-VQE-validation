#!/usr/bin/env python3
"""Compare the local Suzuki commutator heuristic with an amplitude-only baseline."""
import argparse, pandas as pd
from scipy.stats import spearmanr
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("input_csv")
    ap.add_argument("--output",default="source_data_suzuki_diagnostic_comparison_v6.csv")
    a=ap.parse_args(); df=pd.read_csv(a.input_csv); rows=[]
    for case,g in df.groupby("case",sort=False):
        zero=g["comm_relF"]==0; nz=~zero
        rows.append({"case":case,"n_applications":len(g),"n_zero_commutator":int(zero.sum()),
        "max_local_infidelity_when_ccomm_zero":float(g.loc[zero,"local_infidelity"].max()),
        "spearman_theta3_all":float(spearmanr(g["theta_abs"]**3,g["local_infidelity"]).statistic),
        "spearman_theta3_ccomm_all":float(spearmanr(g["theta3_comm"],g["local_infidelity"]).statistic),
        "spearman_theta3_nonzero_ccomm":float(spearmanr(g.loc[nz,"theta_abs"]**3,g.loc[nz,"local_infidelity"]).statistic),
        "spearman_theta3_ccomm_nonzero_ccomm":float(spearmanr(g.loc[nz,"theta3_comm"],g.loc[nz,"local_infidelity"]).statistic)})
    out=pd.DataFrame(rows); out.to_csv(a.output,index=False); print(out.to_string(index=False))
if __name__=="__main__": main()
