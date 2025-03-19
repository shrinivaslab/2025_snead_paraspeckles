import numpy as np
import os
import time
import csv


def input_params(inputfile):
    """
    Reads and stores input parameters from a CSV file into a dictionary.
    Args:
    - inputfile (str): path to the input CSV file
    Returns:
    - input_params (dict): dictionary containing input parameters
    """

    #open file and read parameters
    input_params={}
    with open(inputfile,'r') as csvfile:
        reader=csv.reader(csvfile)
        for row in reader:
            if not row:
                continue
            if row[0][0] !='#':
                key=row[0]
                values=row[1:]

                #catagorize the keys that are floats/integers
                if key in [
                    "Npolymer",
                    "L",
                    "LE",
                    "LM",
                    "vrna",
                    "Lbox",
                    "usemdseed",
                    "seed",
                    "mdseed",
                    "spacing",
                    "Nproteins",
                    "vproteins",
                    "vP1",
                    "vP2",
                    "vP3",
                    "hk",
                    "r0",
                    "r_cut",
                    "Msigma",
                    "Mep",
                    "Esigma",
                    "Eep",
                    "EMsigma",
                    "EMep",
                    "P1ep",
                    "P2ep",
                    "P3ep",
                    "P12ep",
                    "P13ep",
                    "P23ep",
                    "Psigma",
                    "P1Mep",
                    "P2Mep",
                    "P3Mep",
                    "PMsigma",
                    "P1E1ep",
                    "P1E2ep",
                    "P2E1ep",
                    "P2E2ep",
                    "P3E1ep",
                    "P3E2ep",
                    "PEsigma",
                    "initialperiod",
                    "finalperiod",
                    "dt",
                    "firsthalf",
                    "secondhalf",
                ]: 
                    if values:
                        values=sum([float(v) for v in values])
                    else:
                        values=float(values[0])

                #catagorize the keys that are strings
                elif key in ["vcalc","InitialStateFilePath","WalkType","JobInputName","TrajectoryFilePath","DataAnalysis","InterMatrix"]:
                    values=' '.join(values) if values else ''

                #catagorize the keys that are lists
                else: 
                    if isinstance(values,list):
                        values=' '.join(values) if values else ''
                    elif isinstance(values,str):
                        values=' '.join(values) if values else ''
                    else:
                        values=float(values[0])
                input_params[key]=values
    return input_params

def interactionmatrix(inputfile):
    """
    Reads and stores interaction matrix parameters from a CSV file into a dictionary.
    Args:
    - inputfile (str): path to the input CSV file
    Returns:
    - intermatrix (dict): dictionary containing interaction matrix parameters
    """

    #open file and read parameters
    intermatrix={}
    with open(inputfile,'r') as csvfile:
        reader=csv.reader(csvfile)
        for row in reader:
            if not row:
                continue
            if row[0][0] !='#':
                key=row[0]
                values=row[1:]

                #catagorize the keys that are floats/integers
                if key in [
                    "Mep",
                    "Eep",
                    "EMep",
                    "P1ep",
                    "P2ep",
                    "P3ep",
                    "P1Mep",
                    "P2Mep",
                    "P3Mep",
                    "P12ep",
                    "P13ep",
                    "P23ep",
                    "P1E1ep",
                    "P1E2ep",
                    "P2E1ep",
                    "P2E2ep",
                    "P3E1ep",
                    "P3E2ep",
                ]: 
                    if values:
                        values=sum([float(v) for v in values])
                    else:
                        values=float(values[0])
                #catagorize the keys that are lists
                else: 
                    if values:
                        values=sum([float(x) for x in values])
                    else:
                        values=float(values[0])
                intermatrix[key]=values
    return intermatrix


def write_dict_to_csv(dictionary, filename):
    """
    Writes a dictionary to a CSV file.
    Args:
    - dictionary (dict): dictionary to write to CSV
    - filename (str): path to the output CSV file
    """

    #create filepath
    os.makedirs(os.path.dirname(filename),exist_ok=True)

    #write dictionary to CSV
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for key, value in dictionary.items():
            if isinstance(value, list):
                writer.writerow([key] + value)
            else:
                writer.writerow([key, value])

def dict_to_string(dictionary, keys):
    """
    Converts a dictionary to a string with specified keys.
    Args:
    - dictionary (dict): dictionary to convert
    - keys (list): list of keys to include in the string
    Returns:
    - result (str): string representation of the dictionary with specified keys
        """
    return '_'.join(f'{k}_{dictionary[k]}' for k in keys)

def create_filename(typefile,key,values):
    """
    Creates a filename based on the type of file, key, and values.
    Args:
    - typefile (str): type of the file
    - key (dict): dictionary containing the keys
    - values (list): list of values to include in the filename
    Returns:
    - filename (str): generated filename
    """
    if typefile=='RNAandProteins':
        base='RNAandProteins'
    elif typefile=='RNAand3Proteins':
        base='RNAand3Proteins'
    elif typefile=='TranscriptionRNA':
        base='TranscriptionRNA'
    filename=base+'_'+dict_to_string(key,values)
    return filename

def make_filepaths(experiment,subfolder,filename):
    """
    Generates file paths for different types of files based on the experiment, subfolder, and filename.
    Args:
    - experiment (str): name of the experiment
    - subfolder (str): name of the subfolder
    - filename (str): base name of the file
    Returns:
    - filepaths (dict): dictionary containing generated file paths
    """
    date=time.strftime("%Y%m%d")
    types=['InitialStateFilePath','JobInputName','TrajectoryFilePath','DataAnalysis']
    filepaths={}
    for i in types:
        if i=='InitialStateFilePath': #files for the initial state of trajectories
            filepath=f'Data/Experiments/{experiment}/{subfolder}/Inputs/{date}/{filename}_InitialState.gsd'
            filepaths['InitialStateFilePath']=filepath
        elif i=='JobInputName': #files for the input parameters of a trajectory
            filepath=f'Data/Experiments/{experiment}/{subfolder}/JobInputs/{date}/{filename}'
            filepaths['JobInputName']=filepath
        elif i=='TrajectoryFilePath': #files for the trajectory of the simulation
            filepath=f'Data/Experiments/{experiment}/{subfolder}/Outputs/{date}/{filename}_Trajectory.gsd'
            filepaths['TrajectoryFilePath']=filepath
        else: #files for the data analysis of the simulation
            filepath=f'Data/Experiments/{experiment}/{subfolder}/DataAnalysis/{filename}_Data.csv'
            filepaths['DataAnalysis']=filepath
        os.makedirs(os.path.dirname(filepath),exist_ok=True)
    return filepaths

def addmatrix(intermatrix,dictionary,rnarna=False, change=False, rnarnatype=['Mep'], rnarnavalue=[1],testtype=['P23ep'], testvalue=[1]):
        """
        Adds interaction matrix parameters to the dictionary and optionally modifies RNA-RNA or RNA-Protein interaction parameters.
        Args:
        - intermatrix (dict): dictionary containing interaction matrix parameters
        - dictionary (dict): dictionary to update with interaction matrix parameters
        - rnarna (bool): flag to indicate if RNA-RNA interaction parameters should be modified
        - change (bool): flag to indicate if specific interaction parameters should be changed
        - rnarnatype (list): list of RNA-RNA interaction parameter types to modify
        - rnarnavalue (list): list of values for the RNA-RNA interaction parameters
        - testtype (list): list of interaction parameter types to change
        - testvalue (list): list of values for the interaction parameters to change
        Returns:
        - dictionary (dict): updated dictionary with interaction matrix parameters and optional modifications
        """
        for key, val in intermatrix.items():
            dictionary[key]=val
        if rnarna:    
            for i,rnavalue in zip(rnarnatype,rnarnavalue):
                dictionary[i]=rnavalue
        if change:
            for v, testv in zip(testtype,testvalue):
                dictionary[v]=testv
        return dictionary

def create_input(intermatrixfile,version, t, experiment, subfolder, changes, name_key,rnarna=False, change=False, rnarnatype=['Mep'], rnarnavalue=[1],testtype=['P23ep'], testvalue=[1]):
    """
    Creates an input file for a simulation based on the provided parameters and interaction matrix.
    Args:
    - intermatrixfile (str): path to the interaction matrix CSV file
    - version (str): version of the simulation
    - t (str): type of simulation
    - experiment (str): name of the experiment
    - subfolder (str): name of the subfolder
    - changes (dict): dictionary containing parameter changes
    - name_key (dict): dictionary containing keys for the filename
    - rnarna (bool): flag to indicate if RNA-RNA interaction parameters should be modified
    - change (bool): flag to indicate if specific interaction parameters should be changed
    - rnarnatype (list): list of RNA-RNA interaction parameter types to modify
    - rnarnavalue (list): list of values for the RNA-RNA interaction parameters
    - testtype (list): list of interaction parameter types to change
    - testvalue (list): list of values for the interaction parameters to change
    Returns:
    - dictionary (dict): updated dictionary with interaction matrix parameters and optional modifications
    """

    #call base input file
    if version=='RNAand3Proteins_2P_ri':
        inputfile='../../Code/SimulationCode/RNAand3proteins_2P_ri.input'
        print(inputfile)
    else: 
        print('File not found')

    #read input parameters and interaction matrix
    dictionary=input_params(inputfile)
    intermatrix=interactionmatrix(intermatrixfile)
    dictionary=addmatrix(intermatrix,dictionary)
    

    #change parameters is necessary
    for key, value in changes.items():
        dictionary[key]=value

    #create filename
    filename=create_filename(t,dictionary,name_key)

    #create filepaths
    filepaths=make_filepaths(experiment,subfolder,filename)
    for name, path in filepaths.items():
        dictionary[name]=path
        
    #create new input file
    jobinput=dictionary['JobInputName']
    write_dict_to_csv(dictionary, jobinput)
    return dictionary