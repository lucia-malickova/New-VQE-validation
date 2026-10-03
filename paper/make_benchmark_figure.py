#!/usr/bin/env python3
"""Regenerate the benchmark-system overview figure used in the current manuscript."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle

HERE = Path(__file__).resolve().parent

ATOM_COLORS = {'Cu':'#c47b2f','N':'#3d74d1','S':'#e3ba12','H':'#efefef','O':'#e04b4b'}
ATOM_EDGES = {'Cu':'#7b4b1c','N':'#264d8f','S':'#8e7300','H':'#9a9a9a','O':'#8a1f1f'}
SIZE_A = {'Cu':1100,'N':820,'S':980,'H':330,'O':820}

def panel_box(ax):
    rect = FancyBboxPatch((0,0),1,1,
        boxstyle='round,pad=0.008,rounding_size=0.03',
        linewidth=1.5, edgecolor='#bdbdbd', facecolor='#f0f0f0',
        transform=ax.transAxes, clip_on=False, zorder=-100)
    ax.add_patch(rect)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(0,1); ax.set_ylim(0,1)
    for s in ax.spines.values():
        s.set_visible(False)

def draw_atom(ax, x, y, el, size=None, fs=None, label=None, z=3):
    if size is None:
        size = SIZE_A.get(el, 700)
    ax.scatter([x], [y], s=size, c=[ATOM_COLORS[el]],
               edgecolors=[ATOM_EDGES[el]], linewidths=1.4,
               zorder=z, transform=ax.transAxes)
    if label is None:
        label = el
    if fs is None:
        fs = 10 if el != 'H' else 9
    ax.text(x, y, label, transform=ax.transAxes, ha='center', va='center',
            fontsize=fs, weight='bold' if el != 'H' else 'normal', zorder=z+1)

def main():
    fig = plt.figure(figsize=(14.0, 7.6), dpi=150)
    axA = fig.add_axes([0.03,0.06,0.57,0.88])
    axB = fig.add_axes([0.64,0.69,0.33,0.24])
    axC = fig.add_axes([0.64,0.38,0.33,0.24])
    axD = fig.add_axes([0.64,0.07,0.33,0.24])
    for ax in (axA, axB, axC, axD):
        panel_box(ax)

    coords = {
        0:('Cu',0.000,0.000,0.000),
        1:('N',2.000,0.000,0.000),
        2:('H',2.340,0.900,0.200),
        3:('H',2.340,-0.620,0.740),
        4:('H',2.340,-0.280,-0.800),
        5:('N',-1.000,1.732,0.000),
        6:('H',-1.280,2.200,0.800),
        7:('H',-1.550,2.120,-0.620),
        8:('H',-0.120,2.260,0.040),
        9:('S',-0.600,-1.100,1.750),
        10:('H',-1.420,-1.580,2.120),
        11:('S',0.200,-0.500,-2.830),
        12:('H',0.920,-1.180,-3.200),
        13:('H',-0.540,-0.920,-3.280),
    }
    bonds = [(0,1),(1,2),(1,3),(1,4),(0,5),(5,6),(5,7),(5,8),
             (0,9),(9,10),(0,11),(11,12),(11,13)]

    axA.text(0.04,0.965,'(a)  Cu(II) first-shell model', fontsize=20,
             weight='bold', va='top', transform=axA.transAxes)
    axA.text(0.19,0.905,r'$x$--$z$ projection of the archived geometry',
             transform=axA.transAxes, fontsize=10.5, color='#444444', va='center')

    legend_x, yL = 0.58, 0.885
    for k, el in enumerate(['Cu','N','S','H']):
        xx = legend_x + 0.10*k
        draw_atom(axA, xx, yL, el, size=220, fs=8, z=10)
        axA.text(xx+0.022, yL, el, transform=axA.transAxes,
                 va='center', ha='left', fontsize=9)

    proj = np.array([[v[1], v[3]] for v in coords.values()])
    minx,maxx = proj[:,0].min(),proj[:,0].max()
    minz,maxz = proj[:,1].min(),proj[:,1].max()
    left,right,bottom,top = 0.16,0.72,0.24,0.72
    scale = min((right-left)/(maxx-minx), (top-bottom)/(maxz-minz))
    cx,cz = (minx+maxx)/2, (minz+maxz)/2
    cxT,czT = (left+right)/2, (bottom+top)/2
    P = {}
    for i,(el,x,y,z) in coords.items():
        P[i] = (cxT+(x-cx)*scale, czT+(z-cz)*scale, y)

    for i,j in sorted(bonds, key=lambda ij: abs((P[ij[0]][2]+P[ij[1]][2])/2)):
        x1,y1,_ = P[i]; x2,y2,_ = P[j]
        axA.plot([x1,x2],[y1,y2], color='#6f6f6f', lw=2.8,
                 solid_capstyle='round', zorder=1, transform=axA.transAxes)
    for i in sorted(P, key=lambda ii: P[ii][2]):
        el = coords[i][0]
        x,y,_ = P[i]
        draw_atom(axA, x, y, el, size=SIZE_A[el], fs=10 if el!='H' else 9, z=3)

    ox,oy = 0.095,0.81
    axA.annotate('', xy=(ox+0.09,oy), xytext=(ox,oy), xycoords=axA.transAxes,
                 arrowprops=dict(arrowstyle='->', lw=1.35, color='#333333'))
    axA.text(ox+0.10, oy-0.004, 'x', transform=axA.transAxes, fontsize=9, va='center')
    axA.annotate('', xy=(ox,oy-0.09), xytext=(ox,oy), xycoords=axA.transAxes,
                 arrowprops=dict(arrowstyle='->', lw=1.35, color='#333333'))
    axA.text(ox-0.012, oy-0.10, 'z', transform=axA.transAxes,
             fontsize=9, ha='right', va='center')
    cy = oy+0.07
    axA.add_patch(Circle((ox,cy),0.012, transform=axA.transAxes,
                         facecolor='none', edgecolor='#333333', lw=1.2))
    axA.plot([ox],[cy], marker='o', ms=3.2, color='#333333', transform=axA.transAxes)
    axA.text(ox+0.03, cy, 'y', transform=axA.transAxes, fontsize=9, va='center')

    cu=np.array(P[0][:2]); s9=np.array(P[9][:2]); s11=np.array(P[11][:2])
    mid_axial=0.60*cu+0.40*s9
    mid_short=0.58*cu+0.42*s11
    axA.annotate('axial SH$_2$', xy=tuple(mid_axial),
                 xytext=(mid_axial[0]+0.08, mid_axial[1]+0.05),
                 textcoords=axA.transAxes, xycoords=axA.transAxes,
                 arrowprops=dict(arrowstyle='-', color='#777777', lw=1.2),
                 fontsize=11.5, color='#333333')
    axA.annotate('short SH$^{-}$', xy=tuple(mid_short),
                 xytext=(mid_short[0]-0.14, mid_short[1]-0.07),
                 textcoords=axA.transAxes, xycoords=axA.transAxes,
                 arrowprops=dict(arrowstyle='-', color='#777777', lw=1.2),
                 fontsize=11.5, color='#333333')

    bar=2.0*scale
    x0,y0=0.06,0.155
    axA.plot([x0,x0+bar],[y0,y0], color='#333333', lw=2.0, transform=axA.transAxes)
    axA.plot([x0,x0],[y0-0.009,y0+0.009], color='#333333', lw=1.3, transform=axA.transAxes)
    axA.plot([x0+bar,x0+bar],[y0-0.009,y0+0.009], color='#333333', lw=1.3, transform=axA.transAxes)
    axA.text(x0+bar/2,y0+0.018,'2 Å', transform=axA.transAxes,
             ha='center', va='bottom', fontsize=10)

    axA.text(0.04,0.085,'CAS(15e,9o)  •  18 qubits  •  primary molecular benchmark',
             transform=axA.transAxes, fontsize=14.5)
    axA.text(0.04,0.045,'two NH$_3$ donors, short SH$^{-}$ donor, longer axial SH$_2$ donor',
             transform=axA.transAxes, fontsize=11.5, color='#444444')

    def diatomic_panel(ax, title, left_el, right_el, subtext):
        ax.text(0.04,0.92,title, fontsize=18, weight='bold', va='top', transform=ax.transAxes)
        y=0.57; x1=0.32; x2=0.75
        ax.plot([x1,x2],[y,y], color='#6e6e6e', lw=4.0,
                solid_capstyle='round', transform=ax.transAxes, zorder=1)
        draw_atom(ax, x1, y, left_el, size=800, fs=16, z=5)
        draw_atom(ax, x2, y, right_el, size=800, fs=16, z=5)
        ax.text(0.04,0.13,subtext, transform=ax.transAxes, fontsize=15)

    diatomic_panel(axB,'(b)  NO radical','N','O',
                   r'CAS(7e,6o)  •  12 qubits  •  $R_{mathrm{NO}}=1.1508$ Å')
    diatomic_panel(axC,'(c)  OH radical','O','H',
                   r'CAS(7e,5o)  •  10 qubits  •  $R_{mathrm{OH}}=0.970$ Å')

    axD.text(0.04,0.92,'(d)  Open Hubbard chain', fontsize=18,
             weight='bold', va='top', transform=axD.transAxes)
    y=0.56; xs=[0.20,0.43,0.66,0.89]
    for a,b in zip(xs[:-1],xs[1:]):
        axD.plot([a,b],[y,y], color='#6e6e6e', lw=4.0,
                 solid_capstyle='round', transform=axD.transAxes, zorder=1)
    for i,x in enumerate(xs,1):
        axD.scatter([x],[y], s=700, c=['#e8edf4'], edgecolors=['#5b6675'],
                    linewidths=1.5, zorder=3, transform=axD.transAxes)
        axD.text(x,y,str(i), transform=axD.transAxes, ha='center', va='center',
                 fontsize=15, weight='bold', zorder=4)
    axD.text((xs[1]+xs[2])/2, y+0.075, r'$t$', transform=axD.transAxes,
             fontsize=18, ha='center', va='center')
    u_site_x=xs[3]
    axD.plot([u_site_x,u_site_x],[y+0.055,y+0.108], color='#6e6e6e',
             lw=1.8, transform=axD.transAxes)
    axD.text(u_site_x+0.005, y+0.14, r'$U$', transform=axD.transAxes,
             fontsize=18, ha='center', va='center')
    axD.text(0.04,0.13,r'4 sites  •  5 electrons  •  $U/t=1,2,4,8$',
             transform=axD.transAxes, fontsize=15)

    fig.savefig(HERE/'fig_benchmark_systems.pdf', bbox_inches='tight')
    fig.savefig(HERE/'fig_benchmark_systems.png', bbox_inches='tight')
    plt.close(fig)
    print('Wrote', HERE/'fig_benchmark_systems.pdf')

if __name__ == '__main__':
    main()
