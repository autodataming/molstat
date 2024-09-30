from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import numpy as np
import os
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from glob import glob 

mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
model = TSNE(n_components=2, random_state=0, perplexity=30, n_iter=5000)

def fp_from_mol(mol):
    fp = mfpgen.GetFingerprint(mol)
    return fp

def mols2fps(mols_gen):
    fps = []
    totalid = 0
    errid = 0
    for mol in mols_gen:   
        totalid += 1
        try:
            Chem.SanitizeMol(mol)
            fp = fp_from_mol(mol)
            np_array = np.zeros((1, 2048), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp, np_array[0])
            fps.append(np_array[0])
        except Exception as e:
            print("An error occurred:", str(e))
            errid += 1
            smi = Chem.MolToSmiles(mol) if mol else "Invalid Molecule"
            print("error smi:", errid, totalid, smi)    
    return fps 

def tSNE(sdf_files):
    basenames = [i.strip('.sdf') for i in sdf_files]
    fps_legend = {}
    
    legends = []
    all_fps = []
    for sdf_file in sdf_files:
        filename = sdf_file
        basename, file_extension = os.path.splitext(sdf_file)
        legend = basename
        if file_extension == '.sdf':
            mols_gen = Chem.SDMolSupplier(filename, sanitize=False)
            fps_array = mols2fps(mols_gen)
            fps_legend[legend] = np.array(fps_array)
            legends.append(legend)
            all_fps.extend(fps_array)
        else:
            print("%s extension format is not supported" % file_extension)

    all_fps = np.array(all_fps)     
    tsne_results = model.fit_transform(all_fps)

    # 总图
    plt.figure(figsize=(10, 6))
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'cyan', 'magenta']
    idd = 0
    
    for legend in legends:
        color = colors[idd]
        start_idx = sum(len(fps_legend[k]) for k in legends[:idd])
        end_idx = start_idx + len(fps_legend[legend])
        
        x = tsne_results[start_idx:end_idx, 0]
        y = tsne_results[start_idx:end_idx, 1]
        
        plt.scatter(x, y, color=color, label=legend, s=5)
        idd += 1

    plt.title('Compound Space Distribution of Different Libraries')
    plt.xlabel('t-SNE1')
    plt.ylabel('t-SNE2')
    plt.legend()
    plt.grid(False)  # 去掉网格线
    plt.savefig('total_compound_space_distribution.png', transparent=True)  # 保存总图
    # plt.show()

    # 导出每个库的单独图
    for legend in legends:
        plt.figure(figsize=(10, 6))
        color = colors[legends.index(legend)]
        start_idx = sum(len(fps_legend[k]) for k in legends[:legends.index(legend)])
        end_idx = start_idx + len(fps_legend[legend])
        
        x = tsne_results[start_idx:end_idx, 0]
        y = tsne_results[start_idx:end_idx, 1]
        
        plt.scatter(x, y, color=color, s=5)
        plt.title(f'Compound Space Distribution of {legend}')
        plt.xlabel('t-SNE1')
        plt.ylabel('t-SNE2')
        plt.grid(False)  # 去掉网格线
        plt.savefig(f'{legend}_compound_space_distribution.png', transparent=True)  # 保存每个库的图
        plt.close()  # 关闭当前图形以释放内存

if __name__ == '__main__':
    sdf_files = glob("*sdf")
    tSNE(sdf_files)
