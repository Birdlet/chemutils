import os
import numpy as np
import matplotlib.pyplot as plt


def run_gmmpbsa():
    bash_cmd = """
    export MMPBSA=/opt/g_mmpbsa/bin
    export GMX=/opt/gromacs-5.1.2
    export GMXLIB=$GMX/share/gromacs/top
    export PYTHON=python3.7

    export LIG=OUT
    export REC=Protein

    $GMX/bin/gmx -quiet -nobackup grompp -f step5_production.mdp -o step5_1_less.tpr -c step4.3_equilibration.gro -n index.ndx -p topol.top

    # Three Steps Version
    # g_mmpbsa -f ../step5_1_less.xtc -s ../step5_1_less.tpr  -n ../index.ndx -pdie 2  -decomp
    # g_mmpbsa -f ../step5_1_less.xtc -s ../step5_1_less.tpr -i mmpbsa.mdp -n ../index.ndx -nomme -pbsa -decomp
    # g_mmpbsa -f ../step5_1_less.xtc -s ../step5_1_less.tpr -i mmpbsa.mdp -n ../index.ndx -nomme -pbsa -decomp -apol sasa.xvg -apcon sasa_contrib.dat

    # Single Step Version
    printf "%s\n%s\n" $REC $LIG | $MMPBSA/g_mmpbsa -f ../whole_step5_1.xtc -s ../step5_1_less.tpr -n ../full.ndx -i pbsa.mdp -dt 10000 -pdie 2 -pbsa -decomp

    $PYTHON $MMPBSA/MmPbSaStat.py -m energy_MM.xvg -p polar.xvg -a apolar.xvg

    $PYTHON $MMPBSA/MmPbSaDecomp.py -bs -nbs 2000 -m contrib_MM.dat -p contrib_pol.dat -a contrib_apol.dat

    printf "%s\n%s\n" $REC $LIG | $MMPBSA/energy2bfac -s ../step5_1_less.tpr -i energyMapIn.dat
    """
    os.system(bash_cmd)

def plot_dG(full_energy_file="full_energy.dat", figure_save='dG_energy.png', figsize=(14, 9), kcal=True):
    fig = plt.figure(figsize=figsize)

    data = np.loadtxt(full_energy_file, skiprows=2, unpack=True)
    t = data[0]
    if kcal:
        k = 4.8
        unit = ' (kcal/mol)'
    else:
        _k = 1.0
        unit = ' (kJ/mol)'

    # dG bind
    dG_bind = data[-1]/k
    anotation = r'%.2f $\pm$ %.2f%s' % (dG_bind.mean(), dG_bind.std(), unit)

    ax = fig.add_subplot(221)
    fig.subplots_adjust(bottom=0.2)

    # ax.fill_between(t, dG_bind, color="blue", linestyle="-", alpha=0.1)
    ax.plot(t, dG_bind, color="black", linestyle="-", marker='o', linewidth=1.0)
    ax.text(s=anotation, x=t.max()*0.1, y=(dG_bind.max()-dG_bind.min())*0.6+dG_bind.min(), color='black')

    ax.set_xlabel("time $t$ (ps)")
    ax.set_ylabel(r"$\Delta G$  Binding Energy" + unit)
    ax.set_xlim(t.min(), t.max())

    # dG MM
    dG_xx = data[-4]/k
    anotation = r'%.2f $\pm$ %.2f%s' % (dG_xx.mean(), dG_xx.std(), unit)

    ax = fig.add_subplot(222)
    fig.subplots_adjust(bottom=0.2)

    ax.plot(t, dG_xx, color="black", linestyle="-", marker='o', linewidth=1.0)
    ax.text(s=anotation, x=t.max()*0.1, y=(dG_xx.max()-dG_xx.min())*0.6+dG_xx.min(), color='black')

    ax.set_xlabel("time $t$ (ps)")
    ax.set_ylabel(r"$\Delta G_{MM}$" + unit)
    ax.set_xlim(t.min(), t.max())

    # dG polar
    dG_xx = data[-3]/k
    anotation = r'%.2f $\pm$ %.2f%s' % (dG_xx.mean(), dG_xx.std(), unit)

    ax = fig.add_subplot(223)
    fig.subplots_adjust(bottom=0.2)

    ax.plot(t, dG_xx, color="black", linestyle="-", marker='o', linewidth=1.0)
    ax.text(s=anotation, x=t.max()*0.1, y=(dG_xx.max()-dG_xx.min())*0.6+dG_xx.min(), color='black')

    ax.set_xlabel("time $t$ (ps)")
    ax.set_ylabel(r"$\Delta G_{polar}$" + unit)
    ax.set_xlim(t.min(), t.max())

    # dG non-polar
    dG_xx = data[-2]/k
    anotation = r'%.2f $\pm$ %.2f%s' % (dG_xx.mean(), dG_xx.std(), unit)

    ax = fig.add_subplot(224)
    fig.subplots_adjust(bottom=0.2)

    ax.plot(t, dG_xx, color="black", linestyle="-", marker='o', linewidth=1.0)
    ax.text(s=anotation, x=t.max()*0.1, y=(dG_xx.max()-dG_xx.min())*0.6+dG_xx.min(), color='black')

    ax.set_xlabel("time $t$ (ps)")
    ax.set_ylabel(r"$\Delta G_{nonpolar}$" + unit)
    ax.set_xlim(t.min(), t.max())
    
    fig.suptitle(r'MMPBSA $\Delta G$ Binding Energy in MD Simulation')

    fig.savefig(figure_save, dpi=300)


def plot_dG_contrib(contrib_energy_file="final_contrib_energy.dat", figure_save='dG_energy_contrib.png', kcal=True):

    data = np.loadtxt(contrib_energy_file, skiprows=2, unpack=True, dtype=np.str)
    res = data[0][:-1]
    resid = np.array(range(len(res)))
    dG_bind = data[-2][:-1].astype(np.float)
    if kcal:
        k = 4.8
        unit = ' (kcal/mol)'
    else:
        _k = 1.0
        unit = ' (kJ/mol)'

    # dG bind
    dG_bind /= k
    w = len(res)/60*10
    h = 6

    fig = plt.figure(figsize=(w, h))
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.2)

    # ax.fill_between(t, dG_bind, color="blue", linestyle="-", alpha=0.1)
    ax.bar(resid, dG_bind, 0.5, color='black', alpha=0.6)

    ax.set_xlabel("residues")
    ax.set_ylabel(r"$\Delta G$  Binding Energy" + unit)
    ax.set_xlim(max(0, resid.min()-2), resid.max()+2)
    plt.xticks(resid, res, rotation=90)
    
    fig.suptitle(r'MMPBSA $\Delta G$ Binding Energy Contribution')

    fig.savefig(figure_save, dpi=300)
    return fig


if __name__ == '__main__':
    os.environ['QT_QPA_PLATFORM']='offscreen'
    WORKDIR = './charmm-gui-inx/mmpbsa'
    
    CURDIR = os.curdir
    
    os.chdir(WORKDIR)
    plot_dG()
    plot_dG_contrib()
    os.chdir(CURDIR)