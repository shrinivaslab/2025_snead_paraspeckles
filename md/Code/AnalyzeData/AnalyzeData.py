
import hoomd
import gsd.hoomd
import gsd
import freud
import numpy as np
import sys
import glob
sys.path.append('../SimulationCode/')
import CumulativeAverageFile as caf
import DataMangement as dm
import re
import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.patches as patches

class DataAverager:
    """ Class to average data from multiple files
     Attributes:
     - data (dict): dictionary containing data from each file
     - file_count (dict): dictionary containing the number of files each data point was found in
     Methods:
     - add_data(dictionary): adds data from a dictionary to the data attribute
     - compute_average(): computes the average of the data attribute
     """

    def __init__(self):
        self.data = {}
        self.file_count = {}

    def add_data(self, dictionary):
        for key,data in dictionary.items():
            if key not in self.data:
                self.data[key] = np.array(data)
                self.file_count[key] = 1
            else:
                self.data[key] += np.array(data)
                self.file_count[key] += 1

    def compute_average(self):
        averaged_data = {key: value / self.file_count[key] for key, value in self.data.items()}
        return averaged_data


class AverageMatrix:
    """ Class to average matrix data from multiple files
     Attributes:
     - data (list): list containing data from each file as a matrix
     - file_count (int): number of files the data was found in
     Methods:
     - add_data(matrix): adds data from a matrix to the data attribute
     - compute_average(): computes the average of the data attribute
     """
    def __init__(self):
        self.data = []
        self.file_count = 0

    def add_data(self, matrix):
        self.data.append(matrix)
        self.file_count += 1

        
    def compute_average(self):
        averaged_data = np.sum(self.data, axis=0) / self.file_count
        return averaged_data


#select file type
def select_file(folderlocation,datecheck=None,extrachecks=None,pattern=re.compile(r'seed_\d+')):
    """Selects files from a list of folder locations based on the experiment type and date
    Args:
    - folderlocation (string): string with the folder location of data to analyze
    - datecheck (string): string with the date to select
    - extrachecks (list): list of strings with additional checks to select   
    - pattern (re.Pattern): regular expression pattern to match the file name
    Returns:
    - trajdict (dict): dictionary with the selected files"""

    trajdict = {}

    outputpath=os.path.join(folderlocation,'Outputs')
    output_files = glob.glob(os.path.join(outputpath,'**', '*_Trajectory.gsd'), recursive=True)
    jobinputpath=os.path.join(folderlocation,'JobInputs')
    jobinput_files = glob.glob(os.path.join(jobinputpath, '**','*'), recursive=True)
    for file in output_files + jobinput_files:
        dirpath=os.path.dirname(file)
        if datecheck and datecheck not in dirpath:
            continue
        if extrachecks and not all(extra in dirpath for extra in extrachecks):
            continue

        match = pattern.search(file)
        if match:
            key = match.group(0)

            if key not in trajdict:
                trajdict[key] = {}
            if 'Output' in dirpath:
                trajdict[key]['Trajectory'] = file
            elif 'JobInput' in dirpath:
                trajdict[key]['JobInputName'] = file
    return trajdict

def dict_to_matrix_all(ipf,title,folder=False,save=False,annotation=False,oldmodel=False):

    """ Create a heatmap of the interaction matrix 
    Args:
    - ipf (string): string with the path to the input file
    - folder (string): string with the folder location to save the file
    - title (string): string with the title of the plot
    - save (bool): boolean to save the plot
    - annotation (bool): boolean to add annotations to the plot
    """
    # create a dictionary of the data    
    ip=dm.input_params(ipf)
    split = '/Data/'
    data_out = ip['DataAnalysis']
    data_path = data_out.find(split)
    new_path = 'Data/'
    DA_path = os.path.join(new_path, data_out[data_path + len(split):])
    DA_path=DA_path.replace('/DataAnalysis/','/DataAnalysis/InteractionMatrix/')
    os.makedirs(os.path.dirname(DA_path), exist_ok=True)
    
    if oldmodel:
        labels = ["5'", 'M', "3'", 'NONO', 'FUS']
        matrix={}

        matrix["5'"]=np.array([ip['Eep'],ip['EMep'],ip['Eep'],ip['P1E1ep'],ip['P2E1ep']])
        matrix['M']=np.array([ip['EMep'],ip['Mep'],ip['EMep'],ip['P1Mep'],ip['P2Mep']])
        matrix["3'"]=np.array([ip['Eep'],ip['EMep'],ip['Eep'],ip['P1E2ep'],ip['P2E2ep']])
        matrix['NONO']=np.array([ip['P1E1ep'],ip['P1Mep'],ip['P1E2ep'],ip['P1ep'],ip['P12ep']])
        matrix['FUS']=np.array([ip['P2E1ep'],ip['P2Mep'],ip['P2E2ep'],ip['P12ep'],ip['P2ep']])
    else:
        labels = ["5'", 'M', "3'", 'NONO', 'FUS', 'Tdp43']
        matrix={}

        matrix["5'"]=np.array([ip['Eep'],ip['EMep'],ip['Eep'],ip['P1E1ep'],ip['P2E1ep'],ip['P3E1ep']])
        matrix['M']=np.array([ip['EMep'],ip['Mep'],ip['EMep'],ip['P1Mep'],ip['P2Mep'],ip['P3Mep']])
        matrix["3'"]=np.array([ip['Eep'],ip['EMep'],ip['Eep'],ip['P1E2ep'],ip['P2E2ep'],ip['P3E2ep']])
        matrix['NONO']=np.array([ip['P1E1ep'],ip['P1Mep'],ip['P1E2ep'],ip['P1ep'],ip['P12ep'], ip['P13ep']])
        matrix['FUS']=np.array([ip['P2E1ep'],ip['P2Mep'],ip['P2E2ep'],ip['P12ep'],ip['P2ep'], ip['P23ep']])
        matrix['Tdp43']=np.array([ip['P3E1ep'],ip['P3Mep'],ip['P3E2ep'],ip['P13ep'],ip['P23ep'], ip['P3ep']])



    df = pd.DataFrame(matrix, index=labels)
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(df, annot=annotation, cmap='Purples', xticklabels=df.columns, yticklabels=df.index, vmin=0, vmax=3,linecolor='black', linewidth=0.25)
    cbar = ax.collections[0].colorbar
    cbar.set_label('Interaction Strength', fontsize=16)  # Set label size
    cbar.ax.tick_params(labelsize=14)  # Set tick label size
    # Add top labels
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=12)
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

    # Remove tick marks
    ax.tick_params(axis='x', which='both', length=0)
    ax.tick_params(axis='y', which='both', length=0)

    # Increase size of y-axis labels
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=12)

    # Add border around the entire matrix
    # Get the current position of the heatmap
    pos = ax.get_position()

    # Calculate the position and size of the rectangle
    x0, y0, width, height = pos.x0, pos.y0, pos.width, pos.height

    # Create a rectangle
    rect = patches.Rectangle((x0, y0), width, height, linewidth=2, edgecolor='black', facecolor='none', transform=ax.figure.transFigure)

    # Add the rectangle to the plot
    ax.add_patch(rect)

    plt.title(title, fontsize=16)
    if save: 
        if folder:
            DA_path=DA_path.replace(DA_path.split('/')[-1],title)
            os.makedirs(os.path.dirname(DA_path), exist_ok=True)
            plt.savefig(DA_path + title + '_interaction_matrix.png', dpi=300, bbox_inches='tight')
        else:
            plt.savefig('interaction_matrix.png', dpi=300, bbox_inches='tight')
    
    plt.show()


def make_datafile(path,dataname,data,average=True,datechange=None):
    """Save data to a pickle file
    Args:
    - path (string): string with the path to the input file
    - dataname (string): string with the name of the data
    - data (dict): dictionary with the data to save
    - average (bool): boolean to save the data as an average
    - datechange (list): list with the date to change in the path
    """
    input_params = dm.input_params(path)
    split = '/Data/'
    data_out = input_params['DataAnalysis']
    data_path = data_out.find(split)
    new_path = 'Data/'
    DA_path = os.path.join(new_path, data_out[data_path + len(split):])

    if not average:
        filename=dataname+'_dict'
    else:
        filename=dataname+'_avg'
    DA_path=DA_path.replace('_Data.csv',f'_{filename}')
    DA_path=DA_path.replace('/DataAnalysis/','/DataAnalysis/'+dataname+'/')

    if datechange:
        DA_path=DA_path.replace(datechange[0],datechange[1])
    print(DA_path)
    os.makedirs(os.path.dirname(DA_path), exist_ok=True)
    with open(DA_path, 'wb') as f:
        pickle.dump(data, f)


def find_valid_clusters(file,num_files='all'):
    """Find the valid clusters in a file
    Args:
    - file (string): string with the path to the file
    - num_files (int): number of files to read
    Returns:
    - valid_cluster_ids (list): list with the valid cluster ids
    - largest_cluster_id (list): list with the largest cluster ids
    - np.mean(largest_rg) (float): mean of the largest radii of gyration
    - valid_rgs (dict): dictionary with the radii of gyration of the valid clusters
    """
    with gsd.hoomd.open(file['Trajectory'],mode='r') as f:
        RNA_TYPE_ID={0,1,2}
        PROTEIN_TYPE_ID={3,4,5}
        largest_cluster_id=[]
        largest_rg=[]
        if num_files=='all':
            num_files=len(f)
        
        for frame in f[num_files:]:
            cl=freud.cluster.Cluster()
            box = frame.configuration.box
            particles=frame.particles.position
            particle_types = frame.particles.typeid
            total_system= freud.AABBQuery(box, particles)
            cl = freud.cluster.Cluster()
            cl.compute(total_system, neighbors={'r_max': 1.5})
            clp = freud.cluster.ClusterProperties()
            clp.compute(total_system, cl.cluster_idx)
            centers=clp.centers
            rg=clp.radii_of_gyration

            #valid centers
            valid_cluster_ids = []
            for n, center in enumerate(centers):
                cluster_types = particle_types[cl.cluster_idx == n]
                if any(t in RNA_TYPE_ID for t in cluster_types) and any(t in PROTEIN_TYPE_ID for t in cluster_types):
                    valid_cluster_ids.append(n)

            valid_cluster_masses = {n: np.sum(cl.cluster_idx == n) for n in valid_cluster_ids}
            valid_rgs={n:rg[n] for n in valid_cluster_ids}  
            lc=max(valid_cluster_masses, key=valid_cluster_masses.get)
            largest_cluster_id.append(lc)
            largest_rg.append(rg[lc])
            
        return valid_cluster_ids, largest_cluster_id, np.mean(largest_rg), valid_rgs


def avg_PE(trajdict,datechange=None):
    """Average the potential energy of the system
    Args:
    - trajdict (dict): dictionary with the file paths
    - datechange (list): list with the date to change in the path
    """
#potential energy
    PE_dict={}  
    PE_avg=DataAverager()
    for file, paths in trajdict.items():
        if 'Trajectory' not in paths:
            continue
        time_v,PE_v=caf.calc_PE(paths['Trajectory'])
        PE_dict[file]={'time':time_v,'PE':PE_v}
        PE_avg.add_data(PE_dict[file])
    pdata=PE_avg.compute_average()
    make_datafile(paths['JobInputName'],'PE',pdata,average=True,datechange=datechange)
    make_datafile(paths['JobInputName'],'PE',PE_dict,average=False,datechange=datechange)

def avg_cluster(trajdict,datechange=None):
    """Average the number of clusters, the mass of the clusters, and the radii of gyration of the clusters
    Args:
    - trajdict (dict): dictionary with the file paths
    - datechange (list): list with the date to change in the path"""
    cluster_data_dict={}
    cluster_data=DataAverager()
    for file,paths in trajdict.items():
        if 'Trajectory' not in paths:
            continue
        num_clusters,cluster_mass,cluster_rg,times=caf.num_cluster_per_file(paths['Trajectory'])
        cluster_data_dict[file]={'num_clusters_avg':num_clusters,'cluster_mass':cluster_mass,'cluster_rg':cluster_rg,'time':times}
        cluster_data.add_data(cluster_data_dict[file])
    cdata=cluster_data.compute_average()
    make_datafile(paths['JobInputName'],'ClusterData',cdata,average=True,datechange=datechange)
    make_datafile(paths['JobInputName'],'ClusterData',cluster_data_dict,average=False,datechange=datechange)



def avg_MIP(trajdict,num_frames,datechange=None,LC=False):
    """Average the MIP along the x, y, and z axes
    Args:
    - trajdict (dict): dictionary with the file paths
    - num_frames (int): number of frames to read
    - datechange (list): list with the date to change in the path"""

    matrix_dict={}
    matrix_avg=DataAverager()
    sphere_dict={}
    sphere_avg=DataAverager()
    axis_labels=['x','y','z']
    for axis in axis_labels:
        for file,paths in trajdict.items():
            if 'Trajectory' not in paths:
                continue
            if LC:
                valid_clusters,largest_cluster,rg,rgs=find_valid_clusters(paths,-1)
                matricies=caf.MIP(paths['Trajectory'],largest_cluster,num_frames,axis)
            else:
                valid_clusters,largest_cluster,rg,rgs=find_valid_clusters(paths,-1)
                matricies=caf.MIP(paths['Trajectory'],valid_clusters,num_frames,axis)
            matrix_dict[file]=matricies
            matrix_avg.add_data(matrix_dict[file])
        matrixdata=matrix_avg.compute_average()
        if LC:
            name='MIP_LC'
        else:
            name='MIP_AC'
        make_datafile(paths['JobInputName'],name+'_'+axis,matrixdata,average=True,datechange=datechange)
        make_datafile(paths['JobInputName'],name+'_'+axis,matrix_dict,average=False,datechange=datechange)



def den_avg(trajdict,num_frames,datechange=None,LC=False):
    """Average the intensity profile of species in a cluster
    Args:
    - trajdict (dict): dictionary with the file paths
    - num_frames (int): number of frames to read
    - datechange (list): list with the date to change in the path
    
    Note: Neglecting normalized density data"""

    den_dict={}
    den_avg=DataAverager()

    for file,paths in trajdict.items():
        if 'Trajectory' not in paths:
            continue
        if LC:
            valid_clusters,largest_cluster_id,rg,rgs=find_valid_clusters(paths,-num_frames)
            rdf=caf.calc_density(paths['Trajectory'],largest_cluster_id,num_frames)
        else:
            valid_clusters,cluster_id,rg,rgs=find_valid_clusters(paths,-num_frames)
            rdf=caf.calc_density(paths['Trajectory'],valid_clusters,num_frames)
        den_dict[file]=rdf
        den_avg.add_data(den_dict[file])
    dendata=den_avg.compute_average()
    if LC:
        name='Density_LC'  
    else:
        name='Density_AC'

    make_datafile(paths['JobInputName'],name,dendata,average=True,datechange=datechange)
    make_datafile(paths['JobInputName'],name,den_dict,average=False,datechange=datechange)


def den_rna_avg(trajdict,num_frames,datechange=None):
    """Average the density of RNA along the cluster
    Args:
    - trajdict (dict): dictionary with the file paths
    - num_frames (int): number of frames to read
    - datechange (list): list with the date to change in the path"""
    rna_rdf_avg=DataAverager()
    rna_rdf_dict={}
    for file,paths in trajdict.items():
        if 'Trajectory' not in paths:
            continue
        valid_clusters,largest_cluster_id,rg,rgs=find_valid_clusters(paths,-num_frames)
        rna_rdf=caf.calc_density_along_rna(paths['Trajectory'],largest_cluster_id,num_frames)
        rna_rdf_dict[file]=rna_rdf
        rna_rdf_avg.add_data(rna_rdf_dict[file])
    rna_rdfdata=rna_rdf_avg.compute_average()
    make_datafile(paths['JobInputName'],'RNADen_LC',rna_rdfdata,average=True,datechange=datechange)
    make_datafile(paths['JobInputName'],'RNADen_LC',rna_rdf_dict,average=False,datechange=datechange)



def calc_avg(trajdict,calckey,num_files,LC,datechange=None,save=True,oldmodel=False):
    """Calculate the average of the data
    Args:
    - trajdict (dict): dictionary with the file paths
    - calckey (string): string with the type of calculation to perform
    - num_files (int): number of files to read
    - datechange (list): list with the date to change in the path
    - save (bool): boolean to save the data"""
    i=0
    for seed, paths in trajdict.items():
        if i==0:
            dict_to_matrix_all(paths['JobInputName'],'Matrix Inputs',folder=True,save=save,annotation=True,oldmodel=oldmodel)
            i+=1
        else: 
            continue

    for key in calckey:
        if key=='PE':
            avg_PE(trajdict,datechange)
        elif key=='cluster':
            avg_cluster(trajdict,datechange)
        elif key=='MIP':
            avg_MIP(trajdict,num_files,datechange,LC)
        elif key=='Density':
            den_avg(trajdict,num_files,datechange,LC)
        elif key=='RNADensity':
            den_rna_avg(trajdict,num_files,datechange)
        elif key=='matrixinput':
            dict_to_matrix_all(paths['JobInputName'],'','Matrix Inputs',save=save)
        elif key=='all':
            avg_PE(trajdict,datechange)
            avg_cluster(trajdict,datechange)
            avg_MIP(trajdict,num_files,datechange,LC)
            den_avg(trajdict,num_files,datechange,LC)
            den_rna_avg(trajdict,num_files,datechange)
        else:
            print('Invalid calculation key')
            sys.exit(1)


def run_avg(folderlocations,datecheck=None,extrachecks=None,pattern=re.compile(r'seed_\d+'),calckey='all',num_files=100,datechange=None,LC=False,oldmodel=False):
    """Run the average of the data
    Args:  
    - folderlocations (list): list of strings with the folder locations
    - experimenttype (string): string with the experiment type to select
    - datecheck (string): string with the date to select
    - extrachecks (list): list of strings with additional checks to select
    - extratypes (list): list of strings with additional types to select
    - pattern (re.Pattern): regular expression pattern to match the file name
    - calckey (string): string with the type of calculation to perform
    - num_files (int): number of files to read
    - datechange (list): list with the date to change in the path
    - LC (bool): boolean to calculate the data for the largest cluster"""

    trajdict=select_file(folderlocations,datecheck,extrachecks,pattern)
    print(f'Number of trajectories to analyze: {len(trajdict)}')
    calc_avg(trajdict,calckey,num_files,LC,datechange,save=True,oldmodel=oldmodel)