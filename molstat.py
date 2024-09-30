
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np
# import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import os 
import argparse
from tSNE.tSNE import tSNE 


import warnings

# 忽略所有运行时警告
warnings.filterwarnings("ignore", category=RuntimeWarning)
from rdkit import RDLogger

# 关闭RDKit的警告信息
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)  # 设置日志级别为CRITICAL，忽略警告和错误信息



def getDescriptors(mols_gen):
    desc ={}

    tpsas = []
    logps = []
    mws = []
    hbds = []
    hbas = []

    errid=0
    totalid=0
    for mol in mols_gen:   
        totalid+=1
        try:
            # 尝试清理分子，捕获清理过程中的错误
            Chem.SanitizeMol(mol)
            tpsa_m = Descriptors.TPSA(mol)
            logp_m = Descriptors.MolLogP(mol)
            mw_m = Descriptors.MolWt(mol)
            hbd_m = Descriptors.NumHDonors(mol)
            hba_m = Descriptors.NumHAcceptors(mol)
            tpsas.append(tpsa_m)
            logps.append(logp_m)
            mws.append(mw_m)
            hbds.append(hbd_m)
            hbas.append(hba_m)
        except:
            errid+=1
            smi = Chem.MolToSmiles(mol) if mol else "Invalid Molecule"
            print("error smi:",errid,totalid,smi)





    desc['tpsa']=tpsas
    desc['logp']=logps 
    desc['mw']=mws 
    desc['hbd']=hbds
    desc['hda']=hbas    
    return desc 


def mwfigure(mws,outf='mw_stat.png',title="Molecular Weight Distribution",xlabel="Molecular Weight (g/mol)",ylabel="Frequency"):
    bins = [0, 100, 150, 200, 250, 300, 350,400,450, 500, 550, 600, 650, 700, np.inf]
    # Plotting the distribution using Matplotlib
    plt.figure(figsize=(8, 4))
    plt.hist(mws, bins=bins, color='green', edgecolor='black')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(False)
    plt.savefig(outf, dpi=300)  # Save the figure with 300 dpi resolution
    # plt.show()


def molstat(molf):
    '''
    help:
    example: molstat xxx.sdf 
    
    '''
    filename =molf 
    # Get the file name with extension
    file_name_with_ext = os.path.basename(filename)

    # Split the file name and extension
    basename, file_extension = os.path.splitext(file_name_with_ext)
    if file_extension=='.sdf':
        mols_gen = Chem.SDMolSupplier(filename, sanitize=False)
    else:
        print("%s extension format is not supported"%file_extension)
    

    desc={}
    desc = getDescriptors(mols_gen)
    mwoutf = "MW"+basename+".png"
    mwfigure(desc['mw'],mwoutf)


def main():
    parser = argparse.ArgumentParser(description="Calculate molecular descriptors and plot molecular weight distribution.")
    parser.add_argument('filenames', type=str, nargs='+', help='Path(s) to the SDF file(s) containing molecular structures.')
    parser.add_argument('-f', '--function', type=str, choices=['tSNE'], help='Function to apply (currently supports tSNE).')
    
    args = parser.parse_args()
    
    if args.function == 'tSNE':
        print("tSNE function not implemented yet.")  # Placeholder for future implementation
        print("filenames",args.filenames)
        tSNE(args.filenames)
    else:
        molstat(args.filenames[0])

if __name__ == '__main__':
    '''
    '''
    filename = "test.sdf"
    filename = "L3400.sdf"
    # molstat(filename)
    main()
