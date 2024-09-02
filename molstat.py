
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np
# import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import os 
import argparse


def getDescriptors(mols_gen):
    desc ={}

    tpsas = []
    logps = []
    mws = []
    hbds = []
    hbas = []

    for mol in mols_gen:   
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
    plt.show()


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
        mols_gen = Chem.SDMolSupplier(filename)
    else:
        print("%s extension format is not supported"%file_extension)
    

    desc={}
    desc = getDescriptors(mols_gen)
    mwoutf = "MW"+basename+".png"
    mwfigure(desc['mw'],mwoutf)


def main():
    parser = argparse.ArgumentParser(description="Calculate molecular descriptors and plot molecular weight distribution.")
    parser.add_argument('filename', type=str, help='Path to the SDF file containing molecular structures.')
    args = parser.parse_args()
    molstat(args.filename)



if __name__ == '__main__':
    '''
    '''
    filename = "test.sdf"
    molstat(filename)
