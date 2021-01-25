from oddt import interactions
from oddt import toolkit
import numpy as np
import os


interaction_table = ['hbond', 'salt_bridges','halogen_bond', 'pi_stack', 'pi_cation', 'hydrophobic']
color_table = {'hbond': 'gold', 'salt_bridges':'yellow',
                'halogen_bond':'magenta', 'pi_stack':'blue',
                'pi_cation':'cyan', 'hydrophobic':'gray'}
                
def GetInteraction(pdb, lig):
    '''Only Accepte PDB format for protein and MOL2 for Ligand'''
    assert pdb.protein == True
    global interaction_table
    interac = {}
    #Hydrogen Bond:
    hydrogen_bond = interactions.hbonds(pdb, lig, cutoff=4)
    interac['hbond'] = set('%d-%s'%(i,r) for i, r in zip(hydrogen_bond[0]['resnum'], hydrogen_bond[0]['resname']))
    #Salt Bridge
    salt_bridge = interactions.salt_bridges(pdb, lig, cutoff=4)
    interac['salt_bridges'] = set('%d-%s'%(i,r) for i, r in zip(salt_bridge[0]['resnum'], salt_bridge[0]['resname']))
    #Halogen Bond:
    halogen_bond = interactions.halogenbonds(pdb, lig)
    interac['halogen_bond'] = set('%d-%s'%(i,r) for i, r in zip(halogen_bond[0]['resnum'], halogen_bond[0]['resname']))
    #Pi-Pi Stacking:
    pi_stack = interactions.pi_stacking(pdb, lig, cutoff=6)
    interac['pi_stack'] = set('%d-%s'%(i,r) for i, r in zip(pi_stack[0]['resnum'], pi_stack[0]['resname']))
    #Pi-Cation
    pi_cation = interactions.pi_cation(pdb, lig)
    interac['pi_cation'] = set('%d-%s'%(i,r) for i, r in zip(pi_cation[0]['resnum'], pi_cation[0]['resname']))
    #Hyrdophobic Interaion:
    hydrophobic = interactions.hydrophobic_contacts(pdb, lig)
    interac['hydrophobic'] = set('%d-%s'%(i,r) for i, r in zip(hydrophobic[0]['resnum'], hydrophobic[0]['resname']))
    # Pi-metal & Lig-Metal are ommited
    
    return interac


def toIFPs(interactions):
    """Convert interactions to heatmap

    Args:
        interactions: List

    """
    global interaction_table

    interaction_table_k2i = {}
    for i,k in enumerate(interaction_table):
        # here key is interation type, like hbond or salt-bridge
        interaction_table_k2i[k] = i

    n_frames = len(interactions)
    n_dims = len(interaction_table)
    residues = {}

    for frame, interaction in enumerate(interactions):
        for key in interaction_table:
            for res in interaction[key]:
                if not res in residues:
                    # initial a np.zeros
                    residues[res] = np.zeros((n_dims, n_frames), dtype=np.int8)
                    residues[res][interaction_table_k2i[key], frame] = 1
                else:
                    residues[res][interaction_table_k2i[key], frame] = 1

    residues_table = sorted(residues)
    fingerprints = np.array([residues[res] for res in residues_table])
    return np.array(fingerprints), residues_table, np.array(interaction_table)


def plotFPs(fps, residues_table, interaction_table, figsize=(20, 20)):
    """"""
    import matplotlib.pyplot as plt

    nr,ni,nt = fps.shape  # nr: N residues, ni: N interation, nt: N time frame


    fig = plt.figure(figsize=(w, h), dpi=150)
    
    global color_table
    plt.title('Interaction Persistence in MD Simulation', fontdict={'fontsize':36})
    plt.axis('off')

    k = 0
    for i in range(nr):
        for j in range(ni):
            k += 1
            ax = fig.add_subplot(nr*ni, 1, k)
            # fig.subplots_adjust(bottom=0.2)

            ax.bar(range(nt), fps[i][j], 1.0, color=color_table[interaction_table[j]], alpha=0.6)

            if j == 3:  ax.set_ylabel(residues_table[i], size=20, rotation=90)
            ax.yaxis.tick_right()
            ax2 = ax.twinx()
            persistance = np.sum(fps[i][j])/nt

            ax2.set_ylabel('%.2f' % persistance, size=9, rotation=0)
            ax2.spines['right'].set_visible(False)
            ax2.spines['left'].set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.spines['bottom'].set_visible(False)
            ax2.set_xticks([])
            ax2.set_yticks([])

            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)

            if k == 1:
                ax.spines['top'].set_visible(True)
                ax.spines['top'].set_color('gray')
            else:
                ax.spines['top'].set_visible(False)

            if j%ni == ni-1:
                ax.spines['bottom'].set_visible(True)
                ax.spines['bottom'].set_color('gray')
            else:
                ax.spines['bottom'].set_visible(False)

            ax.set_xticks([])
            ax.set_yticks([])

            ax.set_xlim(0, nt)
            ax.set_ylim(0, 1)
            ax.xaxis.label.set_fontsize(0.1)

    # fig.subplots_adjust(bottom=0.1)
    return fig

if __name__ == '__main__':
    import pickle
    import time
    os.environ['QT_QPA_PLATFORM']='offscreen'

    RUNCMD = False
    RUNINT = True
    RUNFIG = True

    REC = 'prot'
    LIG = 'lig'
    LIGID = 'BI4'
    dTime = 1000 # ps  10000/1000=10 frames
    
    # SEP = '-sep'   # seperate pdb?
    SEP = ''         # or not

    WORKDIR = 'charmm-gui-BI4'
    WORKDIR_DUMPS = os.path.join(WORKDIR, 'dumps')
    if not os.path.exists(WORKDIR_DUMPS): os.mkdir(WORKDIR_DUMPS)

    CURDIR = os.curdir
    os.chdir(WORKDIR_DUMPS)

    if RUNCMD:
        dump_frames_cmd = """
        echo Protein |  gmx trjconv -f ../whole_step5_1.xtc -s ../input.gro -dt %(dTime)d -o %(REC)s.pdb %(SEP)s -nobackup -quiet
        echo %(LIGID)s |  gmx trjconv -f ../whole_step5_1.xtc -s ../input.gro -dt %(dTime)d -o %(LIG)s.pdb %(SEP)s -nobackup -quiet
        obabel -ipdb %(LIG)s.pdb -omol2 -O %(LIG)s.mol2
        sed -i 's/HSD/HIS/g' %(REC)s.pdb
        cat %(LIG)s.pdb|grep ENDMDL|wc
        """ % globals()
        os.system(dump_frames_cmd)

    if RUNINT:
        gpdb = toolkit.readfile('pdb', '%s.pdb'%REC)
        glig = toolkit.readfile('mol2', '%s.mol2'%LIG)
        
        md_interactions = []
        pdb = next(gpdb)
        lig = next(glig)
        i = 0
        start = time.time()
        while True:

            try:
                pdb = next(gpdb)
                pdb.protein = True
            except StopIteration:
                pdb = None
            try:
                lig = next(glig)
            except StopIteration:
                lig = None

            if (lig is None) and (pdb is None):
                break
            elif (lig is None) or (pdb is None):
                print('Skip Frame %d for ligand or receptor parsing error' % i)
                # md_interactions.append( GetInteraction(pdb, lig) )
                i += 1
                continue
            else:
                print('Processing Frame %d' % i)
                md_interactions.append( GetInteraction(pdb, lig) )
                i += 1

        dur = int(time.time()-start)
        print('Processing %d frames in %d min %d s' % (i+1, dur//60, dur%60))
        
        fps, residues_table, interaction_table = toIFPs(md_interactions)
        with open('interations.pkl', 'wb') as f:
            pickle.dump((fps, residues_table, interaction_table), f)
    else:
        with open('interations.pkl', 'rb') as f:
            fps, residues_table, interaction_table = pickle.load(f)
    
    if RUNFIG:
        WIDTH = 20

        nr,ni,nt = fps.shape  # nr: N residues, ni: N interation, nt: N time frame
        h = int(nr*ni/nt * WIDTH * 2)
        w = WIDTH
        if abs(h-w)>w: h = 1.2*w

        fig = plotFPs(fps, residues_table, interaction_table, (w, h))
        fig.savefig('interaction_fps.png', dpi=150)

    os.chdir(CURDIR)