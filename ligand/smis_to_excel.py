import os
from io import BytesIO
import sys

try:
    import xlsxwriter
except ImportError:
    print('ImportError: xlsxwriter is not installed or import failed', file=sys.stderr)
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
except ImportError:
    print('ImportError: rdkit is not installed or import failed', file=sys.stderr)


class XLSXConverter():
    def __init__(self, xlsfile):
        self.work_book = xlsxwriter.Workbook(xlsfile)

    def add_sheet(self, csv, sheet_name=None, header=True, smi_idx=0, sep=',', img_size=(300,300)):
        if not sheet_name:
            sheet_name = os.path.basename(csv)
        sheet = self.work_book.add_worksheet(sheet_name)

        if not os.path.exists(csv):
            print('CSV file "%s" is not found!' % csv)

        with open(csv) as f:
            data = []
            for l in f.readlines():
                data.append( l.split(sep) )
        if header:
            title = data[0]
        else:
            title = ['COL%d'%(i+1) for i in range(len(data[0]))]
        ncol = len( title )
        data = [d for d in data if len(d)==ncol]

        # insert image at col 1st
        sheet.write(0, 0, 'Molecule')
        for i, t in enumerate(title):
            sheet.write(0, i+1, t)
        sheet.set_column("A:A", width=img_size[1]/6)
        sheet.set_default_row(img_size[0]/1.2)
        sheet.set_row(0, 25)
        
        for i, d in enumerate(data):
            # print(i)
            smi = d[smi_idx]
            # print(smi)
            mol = Chem.MolFromSmiles(smi)

            # insert image at col 1st
            if mol:
                try:
                    from rdkit.Chem import rdDepictor
                    rdDepictor.SetPreferCoordGen(True)
                except:
                    continue
                image_data = Draw.MolToImage(mol, size=img_size)
                image_bytes = BytesIO()
                image_data.save(image_bytes, 'png')
                sheet.insert_image(i+1, 0, "%s_smiles_%s.png" % (csv, str(i+1)),
                    {'image_data': image_bytes, "positioning": 1})
            else:
                sheet.write(i+1, 0, 'Convert SMILES to Image Failed')

            for j, c in enumerate(d):
                sheet.write(i+1, j+1, c)
    
    def close(self):
        self.work_book.close()

    def __del__(self):
        self.work_book.close()

def check_resolution(restr):
    import re
    if re.match('^[0-9]+x[0-9]+$', restr):
        return tuple(int(x) for x in restr.split('x'))
    else:
        return None


if __name__ == '__main__':
    import sys
    import argparse


    parser = argparse.ArgumentParser("Draw Molecule Figures in Excel")
    parser.add_argument("-i", "--input", nargs='+', type=str, help="Input CSV file")
    parser.add_argument("-o", "--output", default='out.xls', help="Output XLS file")
    parser.add_argument("-s", "--sep", default=",", help="CSV file seperator")
    parser.add_argument("-x", "--index", default=0, type=int, help="SMILES index column, index start from 0")
    parser.add_argument("--noheader", action="store_false", help="No header line in CSV file")
    parser.add_argument("--figsize", default='300x300', help="Molecule figure size")

    args = parser.parse_args()
    
    if check_resolution(args.figsize) is None:
        figsize=(300,300)
        print('Resolution Setting Failed ( user set: "%s" ), use 300x300 instead' % args.figsize, file=sys.stderr)
    else:
        figsize = check_resolution(args.figsize)
    
    if not args.input:
        print('Input CSVs are not provided, Failed!', file=sys.stderr)
        exit(1)
    
    xls = XLSXConverter(xlsfile=args.output)
    snames = set()  # xlsxwriter.exceptions.DuplicateWorksheetName
    for j, csv in enumerate(args.input):
        sname = os.path.basename(csv)
        if sname in snames:
            sname += '_%d'%(j+1)
        snames.add(sname)
        if not os.path.exists(csv):
            print('Input CSV file %s is not existed, skip' % csv)
            continue
        print('Processing No.%d %s...' % ((j+1), csv))
        xls.add_sheet(csv=csv, sheet_name=sname,smi_idx=args.index, header=args.noheader, sep=args.sep, img_size=figsize)
    xls.close()
    print('Convert Finished! Save %s' % args.output)

