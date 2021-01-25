import requests
import gzip
import os
import time
import re

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

def fetch_pdb_file(pdb_code, file_dest=None):
    if not file_dest:
        file_dest = file_dest + '.pdb'
    print("https://files.rcsb.org/download/%s.pdb.gz" % pdb_code.lower())
    response = requests.get("https://files.rcsb.org/download/%s.pdb.gz" % pdb_code.lower())
    if response.status_code == 404:
        raise ConnectionError("This PDB code: %s does not exist in RCSB PDB, please check this manually by user." % pdb_code)

    pdb_data = gzip.decompress(response.content).decode()
    with open(file_dest, 'w') as f:
        f.write(pdb_data)
    

def fetch_pdb_by_uniprot(acc, runtime='.'):
    uniprot_data = requests.get('https://www.uniprot.org/uniprot/%s.txt' % acc.upper()).text

    seq_flag = 0
    seq = ''
    for l in uniprot_data.split('\n'):
        if l.startswith('DR   PDB;'):
            pdb_row = [c.strip() for c in l.split(';')]
            pdb_code = pdb_row[1]
            file_dest = os.path.join(runtime, pdb_code+'.pdb')
            print(file_dest)
            fetch_pdb_file(pdb_code, file_dest)
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

def fetch_pdb_by_code_file(codefile, runtime='.'):
    with open(codefile) as f:
        for l in f.readlines():
            pdb_code = l.strip().upper()
            if not is_pdb_code(pdb_code): continue
            file_dest = os.path.join(runtime, pdb_code+'.pdb')
            print(file_dest)
            fetch_pdb_file(pdb_code, file_dest)
            time.sleep(2)


def fetch_pdb_by_code(code, runtime='.'):
    file_dest = os.path.join(runtime, code+'.pdb')
    print(file_dest)
    fetch_pdb_file(code, file_dest)
    time.sleep(2)



if __name__ == '__main__':
    import sys
    try:
        acc_code = sys.argv[1]
    except IndexError:
        print('Please provide PDB code, UniProt Accession, or PDB code file is not provided')
        exit()
    try:
        runtime = sys.argv[2]
    except IndexError:
        runtime = '.'
        print('Using default runtime directory "."')
    if is_uniprot_acc(acc_code):
        fetch_pdb_by_uniprot(acc=acc_code, runtime=runtime)
    elif is_pdb_code(acc_code):
        fetch_pdb_by_code(code=acc_code, runtime=runtime)
    elif os.path.isfile(acc_code):
        fetch_pdb_by_code_file(codefile=acc_code, runtime=runtime)
    else:
        print('Please provide PDB code, UniProt Accession, or PDB code file')
