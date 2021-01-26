import requests
import gzip
import os
import time
import re
import json

SOLVENTS = ["EOH", "MOH", "GOL", "1GP", "2DP", "3PH", "6A8", "DGA",
            "DGD", "DGG", "DR9", "DVG", "G3P", "HGP", "HGX", "IGP",
            "INB", "L1P", "L2P", "L3P", "L4P", "LHG", "LI1", "LIO",
            "LPC", "PGM", "SGL", "SGM", "SQD", "TGL", "DMS", "EDO",
            "12P", "15P", "1PE", "2PE", "CE9", "CP4", "DIO", "P4C",
            "P6G", "PG4", "PGE", "VNY", "PEG", "TRS", "IPA", "TBU",
            "ACT", "EEE", "ACY", "BME", "MBN", 
            "SO4", "SO3", "SO2",
            "UNK"]

def is_uniprot_acc(uniprot_acc):
    uniprot_acc = uniprot_acc.upper()
    uniprot_acc_pattr = re.compile('^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$')
    if re.match(uniprot_acc_pattr, uniprot_acc):
        return True
    else:
        return False

def is_pdb_code(pdb_code):
    pdb_code = pdb_code.upper()
    pdb_code_pttr = re.compile('^[1-9][A-Z0-9][A-Z0-9][A-Z0-9]$')
    if re.match(pdb_code_pttr, pdb_code):
        return True
    else:
        return False


def fetch_pdb_by_uniprot(acc, runtime='.'):
    uniprot_data = requests.get('https://www.uniprot.org/uniprot/%s.txt' % acc.upper()).text

    seq_flag = 0
    seq = ''
    pdb_codes = []
    for l in uniprot_data.split('\n'):
        if l.startswith('DR   PDB;'):
            pdb_row = [c.strip() for c in l.split(';')]
            pdb_code = pdb_row[1]
            file_dest = os.path.join(runtime, pdb_code+'.pdb')
            pdb_codes.append(pdb_code)
            time.sleep(2)

        elif l.startswith('SQ'):
            seq_flag = 1
            continue

        elif l.startswith('//'):
            seq_flag = 0
            break

        if seq_flag:
            for c in l:
                if c == ' ': continue
                seq += c
    return pdb_codes

def fetch_pdb_by_code_file(codefile, runtime='.'):
    pdb_codes = []
    with open(codefile) as f:
        for l in f.readlines():
            pdb_code = l.strip().upper()
            if not is_pdb_code(pdb_code): continue
            pdb_codes.append(pdb_code)
    return pdb_codes


# Identify Pocket
def fetch_pocket_pdbekb(pdb_code):
    """https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/:pdbId

    ...
    {
        'chem_comp_name': 'PHOSPHOAMINOPHOSPHONIC ACID-ADENYLATE ESTER',
        'entity_id': 2,
        'residue_number': 1,
        'author_residue_number': 1477,
        'chain_id': 'A',
        'alternate_conformers': 0,
        'author_insertion_code': '',
        'chem_comp_id': 'ANP',
        'struct_asym_id': 'C'},
        ....
    ...

    """
    global SOLVENTS
    pdb_code = pdb_code.lower()
    response = requests.get("https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/%s" % pdb_code)
    if response.status_code == 404:
        print("Fetch Pocket Data Error, because PDBe didn't record any small molecule ligands for this PDB: %s" % pdb_code)
        return None
    elif response.status_code != 200:
        print("Fetch Pocket Data Error, Status Code: %d" % response.status_code)
        return None

    js = json.loads(response.content)
    
    # for pdbe-kb api, deprecated!
    # ligands = {} # ccd, site_residues, pdbe_recommend_chain_id
    # ions = {}
    # for site_for_ligand in js[pdb_code]["data"]:
    #     ccd, site_residues, pdbe_recommend_chain_id = _extract_pocket(site_for_ligand)
    #     if ccd in SOLVENTS:
    #         continue
    #     if len(ccd) < 3:
    #         ions[ccd] = ( site_residues, pdbe_recommend_chain_id )
    #         continue

    #     ligands[ccd] = ( site_residues, pdbe_recommend_chain_id )

    # for pdbe api
    ligands = {} # ccd, chain_id
    ions = {}
    for site_for_ligand in js[pdb_code]:
        ccd, chain_id = site_for_ligand['chem_comp_id'], site_for_ligand['chain_id']
        if ccd in SOLVENTS:
            continue
        if len(ccd) < 3:
            ions[ccd] = ( None, chain_id )  # site_residue, chain_id, where site_residue is None for PDBe API
            continue
        if ccd in ligands:
            ligands[ccd][1].append(chain_id)
        else:
            ligands[ccd] = ( None, [chain_id] )

    return ligands

def get_smiles_by_ccd(ccd):
    with open('database/ligandExpo.smi') as f:
        for l in f:
            c = l.split('\t')
            if len(c) > 1:
                # SMILES CCD NAME
                if c[1] == ccd:
                    return c[0]
    # no CCD found in database
    # use online server
    try:
        import xml.dom.minidom as minidom
        xmlstr = requests.get('https://files.rcsb.org/ligands/download/%s.xml' % ccd.upper()).text
        tree = minidom.parseString(xmlstr)
        root = tree.documentElement
        for elem in root.getElementsByTagName('PDBx:pdbx_chem_comp_descriptor'):
            if elem.getAttribute('type') == 'SMILES_CANONICAL' and elem.getAttribute('program').startswith('OpenEye'):
                smiles = elem.childNodes[1].childNodes[0].data
                return smiles
    except Exception as e:
        print('Fetch ligand smiles from online server failed')
        print(e)
        return None


if __name__ == '__main__':
    import sys
    try:
        acc_code = sys.argv[1]
    except IndexError:
        print('Please provide PDB code, UniProt Accession, or PDB code file is not provided')
        exit()

    if is_uniprot_acc(acc_code):
        pdb_codes = fetch_pdb_by_uniprot(acc=acc_code)
    elif is_pdb_code(acc_code):
        pdb_codes = [acc_code,]
    elif os.path.isfile(acc_code):
        pdb_codes = fetch_pdb_by_code_file(codefile=acc_code)
    else:
        print('Please provide PDB code, UniProt Accession, or PDB code file')

    smiles = []
    for pdb_code in pdb_codes:
        ligands = fetch_pocket_pdbekb(pdb_code)
        if not ligands:
            print('No Ligand CCD for PDB: %s, CONTINUE' % pdb_code)
            continue
        for ccd in ligands:
            smiles.append((pdb_code, get_smiles_by_ccd(ccd)))

    for pdb_smi in smiles:
        print('%s %s' % pdb_smi)
