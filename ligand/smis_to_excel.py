import os
import xlsxwriter
from io import BytesIO
# from urllib.request import urlopen
# import requests
# import json
import traceback

from rdkit import Chem
from rdkit.Chem import Draw


class XLSXConverter():
    def __init__(self, xlsfile):
        self.work_book = xlsxwriter.Workbook(xlsfile)

    def add_sheet(self, csv, sheet_name=None, smi_idx=0, img_size=(300,300)):
        if not sheet_name:
            sheet_name = os.path.basename(csv)
        sheet = self.work_book.add_worksheet(sheet_name)

        if not os.path.exists(csv):
            print('CSV file "%s" is not found!' % csv)

        with open(csv) as f:
            data = []
            for l in f.readlines():
                data.append( l.split(',') )
        title = data[0]
        ncol = len(title)
        data = [d for d in data[1:] if len(d)==ncol]

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


if __name__ == '__main__':
    import sys
    out_xls = sys.argv[1]
    csvs = sys.argv[2:]
    w = XLSXConverter(xlsfile=out_xls)
    for csv in csvs:
        w.add_sheet(csv=csv)
    w.close()