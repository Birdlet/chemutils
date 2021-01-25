from rdkit import Chem

def isNinAromatic(smi):
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        print('RDkit parsing failed')
        return False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            if atom.GetIsAromatic(): return True
            
    return False
    
    
def isNinAromaticWithNeighborH(smi):
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        print('RDkit parsing failed')
        return False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetIsAromatic():
            for atom in atom.GetNeighbors():
                if (atom.GetNumExplicitHs() + atom.GetNumImplicitHs()) > 0:
                    return True
                    
    return False
    
    
if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print('Usage: python aromatic.py SMILES_STRING')
        exit(0)
    smi = sys.argv[1]
    print('isAromatic: ',  isNinAromatic(smi), '  isAromaticNeighborH: ', isNinAromaticWithNeighborH(smi))