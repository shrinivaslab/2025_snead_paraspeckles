
import numpy as np
import os
import pickle
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

#create color list
colorfolder='../../Code/AnalyzeData/ColorforSpecies/'
mip_color_list={}
for colorfile in os.listdir(colorfolder):
    color_data=pickle.load(open(colorfolder+'/'+colorfile,'rb'))
    name=colorfile.split('.')[0]
    mip_color_list[name]=color_data

def select_avg_file(base_folder,subfolders,datecheck=None):
    """Select the average file from the subfolders
    Args:
        base_folder (str): The base folder where the subfolders are located
        subfolders (str): The subfolder to search for the average file
        datecheck (str): The date to search for in the subfolder
    Returns:
        str: The path to the average file
    """
    folder_path = os.path.join(base_folder, subfolders)
    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        for dirpath, _, filenames in os.walk(folder_path):
            if datecheck is not None:
                if datecheck in dirpath:
                    for filename in filenames:
                        if filename.endswith('_avg'):
                            return os.path.join(dirpath, filename), dirpath
            else:
                for filename in filenames:
                    if filename.endswith('_avg'):
                        return os.path.join(dirpath, filename), dirpath
    return None,None

def format_axes(ax, title=None, xlabel=None, ylabel=None, grid=True, fontsize=12):
    """Format the axes of the plot"""
    if title:
        ax.set_title(title, fontsize=fontsize)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=fontsize)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=fontsize)
    if grid:
        ax.grid(True)

def plot_cluster_data(datafile,dirpath,name,save=False):
    """Plot the cluster data that includes the number of clusters, average cluster mass, and average cluster radius of gyration
    Args:
        datafile (str): The path to the data file
        dirpath (str): The directory path to save the plot
        name (str): The name of the plot
        save (bool): Whether to save the plot
    Returns:
        None
    """
    #load the data
    data=pickle.load(open(datafile,'rb'))

    #prepare the plot
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(name, fontsize=16)

    #plot the average number of clusters
    ax[0].plot(data['time'], data['num_clusters_avg'])
    format_axes(ax[0], title='Number of Clusters', xlabel='Time', ylabel='Number of Clusters', fontsize=14)

    #plot the average cluster mass
    ax[1].plot(data['time'], data['cluster_mass'])
    format_axes(ax[1], title='Average Cluster Mass', xlabel='Time', ylabel='Cluster Mass', fontsize=14)

    #plot the average cluster radius of gyration
    ax[2].plot(data['time'], data['cluster_rg'])
    format_axes(ax[2], title='Average Cluster Radius of Gyration', xlabel='Time', ylabel='Cluster Rg', fontsize=14)

    #save the plot
    if save:
        path=dirpath+'/Plots/'
        os.makedirs(path,exist_ok=True)
        plt.savefig(path + name+'_ClusterAnalysis.png')
    plt.show()


def plot_PE(datafile,dirpath,name,save=False):
    """Plot the potential energy of the system
    Args:
        datafile (str): The path to the data file
        dirpath (str): The directory path to save the plot
        name (str): The name of the plot
        save (bool): Whether to save the plot
    Returns:
        None
    """

    #load data
    data=pickle.load(open(datafile,'rb'))

    #prepare the plot
    fig,ax=plt.subplots(1,1,figsize=(10,10))
    fig.suptitle(name, fontsize=16)

    #plot the potential energy
    ax.plot(data['time'], data['PE'])
    format_axes(ax, title='Potential Energy', xlabel='Time', ylabel='Potential Energy', fontsize=14)

    #save the plot
    if save:
        path=dirpath+'/Plots/'
        os.makedirs(path,exist_ok=True)
        plt.savefig(path+ name+'_PotentialEnergy.png')
    plt.show()


def plot_MIP(datafile,dirpath,name,save=False,oldmodel=False):
    """Plot the Maximum Intensity Projection based on the axis
    Args:
        
        datafile (str): The path to the data file
        dirpath (str): The directory path to save the plot
        name (str): The name of the plot
        save (bool): Whether to save the plot
        oldmodel (bool): Whether the model is the old model
    Returns:
        None
    """
    #determine the axis
    if datafile.endswith('MIP_x_avg') or datafile.endswith('MIP_LC_x_avg'):
        axis='YZ Plane' #title of the plot

    elif datafile.endswith('MIP_y_avg') or datafile.endswith('MIP_LC_y_avg'):
        axis='XZ Plane'
    else:
        axis='XY Plane'

    #determine the color and name based on the species
    color_dict = {'M': mip_color_list['NEAT1-M'], 'E1': mip_color_list['NEAT1-E1'], 'E2': mip_color_list['NEAT1-E2'], 'P1':  mip_color_list['NONO'], 'P2': mip_color_list['FUS'], 'P3': mip_color_list['TDP43']}
    name_dict={'total':'Entire Cluster','M':'M','E1':"5'",'E2':"3'",'P1':'NONO','P2':'FUS','P3':'Tdp43'}

    #load the data
    data=pickle.load(open(datafile,'rb'))

    #plot the MIP based on old or new model
    if oldmodel:
        fig,ax=plt.subplots(1,len(data)-2,figsize=(20,4))
        n=0
        for sp,mip in data.items():
            if sp=='total' or sp=='P3':
                continue
            im=ax[n].imshow(mip.T/np.max(mip.T), cmap=color_dict[sp], vmin=0, origin='lower',interpolation='nearest') 
            ax[n].set_xticks([])    
            ax[n].set_yticks([])
            ax[n].set_title(name_dict[sp], fontsize=14)
            divider = make_axes_locatable(ax[n])
            cax = divider.append_axes("bottom", size="5%", pad=0.08)
            cax.xaxis.set_ticks_position('bottom')
            cax.xaxis.set_label_position('bottom') 
            n+=1
        fig.suptitle(axis, fontsize=14)
        fig.supxlabel('Normalized Volume Density', fontsize=14)
        fig.tight_layout()

        if save:
            path=dirpath+'/Plots/'
            os.makedirs(path,exist_ok=True)
            plt.savefig(path+'_'+name+'_'+axis+'_MIP.png',dpi=600)
        plt.show()

    else:
        fig,ax=plt.subplots(1,len(data)-1,figsize=(20,4))
        n=0
        for sp,mip in data.items():
            if sp=='total':
                continue
            im=ax[n].imshow(mip.T/np.max(mip.T), cmap=color_dict[sp], vmin=0, origin='lower',interpolation='nearest') 
            ax[n].set_xticks([])    
            ax[n].set_yticks([])
            ax[n].set_title(name_dict[sp], fontsize=14)
            divider = make_axes_locatable(ax[n])
            cax = divider.append_axes("bottom", size="5%", pad=0.08)
            cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
            #cbar.set_label(f'NVD of {name_dict[sp]}', fontsize=14)
            cax.xaxis.set_ticks_position('bottom')
            cax.xaxis.set_label_position('bottom') 
            n+=1
        fig.suptitle(axis, fontsize=14)
        fig.supxlabel('Normalized Volume Density', fontsize=14)
        fig.tight_layout()
     
        #save the plot
        if save:
            path=dirpath+'/Plots/'
            os.makedirs(path,exist_ok=True)
            plt.savefig(path+'_'+name+'_'+axis+'_MIP.png',dpi=600)
        plt.show()


def plot_den_by_species(datafile,dirpath,name,save=False,oldmodel=False,LC=False):
    """Plot the density of each species
    Args:
        datafile (str): The path to the data file
        dirpath (str): The directory path to save the plot
        name (str): The name of the plot
        save (bool): Whether to save the plot
        oldmodel (bool): Whether the model is the old model
        LC (bool): Whether to plot the largest cluster data
    Returns:
        None"""
    
    #load the data
    data=pickle.load(open(datafile,'rb'))

    #determine the color and name based on the species
    color_dict={'M_rdf':'deeppink','E1_rdf':'yellowgreen','E2_rdf':'darkgreen','P1_rdf':'#8C510A','P2_rdf':'#5E918B','P3_rdf':'#805D94'}
    name_dict={'M_rdf':'M','E1_rdf':"5'",'E2_rdf':"3'",'P1_rdf':'NONO','P2_rdf':'FUS','P3_rdf':'TDP43'}

    #prepare the plot
    figure, axes = plt.subplots(2, 1, figsize=(4, 6))

    #plot the density of each species
    for key in data.keys():
        if key.endswith('rdf'):
            if key in {'M_rdf', 'E1_rdf', 'E2_rdf'}:
                axes[0].plot(data[key.split('_')[0]+'_bin'][-1],data[key][-1], label=name_dict[key],color=color_dict[key],linewidth=2)
                axes[0].set_title('RNA Profile')
                axes[0].legend()
                axes[0].set_xlim(0,1)
            else:
                if oldmodel:
                    if key=='P3_rdf':
                        continue
                    axes[1].plot(data[key.split('_')[0]+'_bin'][-1],data[key][-1], label=name_dict[key],color=color_dict[key],linewidth=2)
                    axes[1].set_title('Protein Profile')
                    axes[1].legend()
                    axes[1].set_xlim(0,1)
                else:
                    axes[1].plot(data[key.split('_')[0]+'_bin'][-1],data[key][-1], label=name_dict[key],color=color_dict[key],linewidth=2)
                    axes[1].set_title('Protein Profile')
                    axes[1].legend()
                    axes[1].set_xlim(0,1)

    #save the plot
    figure.tight_layout(rect=[0.05, 0.05, 0.90, 0.90])
    figure.supxlabel('Normalized Distance from Center')
    figure.supylabel('Number Density')
    figure.subplots_adjust(top=0.85)
    if save:
        path=dirpath+'/Plots/'
        print(f'Saving to {path}')
        os.makedirs(path,exist_ok=True)
        if LC:
            plt.savefig(path+'_'+name+'_DensityperSpeciesLC.png',dpi=600, bbox_inches='tight')
        else:
            plt.savefig(path+'_'+name+'_DensityperSpecies.png',dpi=600, bbox_inches='tight')
    plt.show()


def plot_rna_density(datafile,dirpath,name,save=False):
    """Plotting density of each third of RNA
    Args:
        datafile (str): The path to the data file
        dirpath (str): The directory path to save the plot
        name (str): The name of the plot
        save (bool): Whether to save the plot
    Returns:
        None"""
    
    #define the names of the thirds
    name_dict={'first_rdf':'1-15','middle_rdf':'16-30','last_rdf':'31-45'}

    #load the data
    data=pickle.load(open(datafile,'rb'))

    #prepare the plot
    figure, ax = plt.subplots(1, 1, figsize=(5, 5))

    #plot the density of each third
    for key in data.keys():
        if key.endswith('rdf'):
            x_name = key.replace('rdf', 'bin')
            ax.plot(data[x_name][-1],data[key][-1],label=name_dict[key],linewidth=2)
            ax.set_title(f'RNA Density for {name}')
    ax.legend(title='Sections of NEAT1 Monomers')
    ax.set_title('RNA Density')

    #save the plot
    figure.tight_layout(rect=[0.05, 0.05, 0.90, 0.90])
    figure.supxlabel('Normalized Distance from Center')
    figure.supylabel('Number Density')
    figure.subplots_adjust(top=0.85)
    if save:
        path=dirpath+'/Plots/'
        os.makedirs(path,exist_ok=True)
        plt.savefig(path+'_'+name+'_RNA_Thirds_Density.png',dpi=600, bbox_inches='tight')
    plt.show()



def plot_data(folder,name,plot_LC=True,save=False,oldmodel=False,datecheck=None):
    """
    Plot data from the specified folder based on the given datecheck and name.
    Args:
        folder (str): The folder containing the data.
        datecheck (str): The date to filter the subfolders.
        name (str): The name of the plot.
        plot_LC (bool): Whether to plot the largest cluster data.
        save (bool): Whether to save the plots.
        oldmodel (bool): Whether the model is the old model.
    Returns:
        None
    """

    #iterate through the subfolders
    for subfolders in os.listdir(folder):
        if os.path.isdir(os.path.join(folder,subfolders)):
            avg_file,dirpath=select_avg_file(folder,subfolders,datecheck)
            if avg_file is not None:
                #check if the data is for the largest cluster or not
                if plot_LC:
                    if "AC" in avg_file:
                        continue 
                else:
                    if "LC" in avg_file:
                        continue
                
                #plot the data based on the file name
                if os.path.exists(avg_file):
                    if avg_file.endswith('ClusterData_avg'):
                        plot_cluster_data(avg_file,dirpath,name,save)
                    elif avg_file.endswith('PE_avg'):
                        plot_PE(avg_file,dirpath,name,save)
                    elif avg_file.endswith('MIP_AC_x_avg') or avg_file.endswith('MIP_AC_y_avg') or avg_file.endswith('MIP_AC_z_avg'):
                        plot_MIP(avg_file,dirpath,name,save,oldmodel)
                    elif avg_file.endswith('MIP_LC_x_avg') or avg_file.endswith('MIP_LC_y_avg') or avg_file.endswith('MIP_LC_z_avg'):
                        plot_MIP(avg_file,dirpath,name,save,oldmodel)
                    elif avg_file.endswith('Density_AC_avg'):
                        plot_den_by_species(avg_file,dirpath,name,save,oldmodel)
                    elif avg_file.endswith('Density_LC_avg'):
                        plot_den_by_species(avg_file,dirpath,name,save,oldmodel,LC=True)
                    elif avg_file.endswith('RNADen_LC_avg'):
                        plot_rna_density(avg_file,dirpath,name,save)
                else:
                    print(f"Data file not found for {avg_file}")
            else:
                print(f"Average file not found for {subfolders}")


def plot_den_alltrajs_density(trajs, avg, oldmodel=False, std=True):
    """Plot intensity by species with optional standard deviation.

    Args:
    - trajs (str): Path to the trajectory data file (pickle format).
    - avg (str): Path to the average data file (pickle format).
    - oldmodel (bool): Whether to exclude certain species (e.g., P3) for old models.
    - std (bool): Whether to plot standard deviation (default: True).
    """
    # Load data
    data = pickle.load(open(trajs, 'rb'))
    avgdata = pickle.load(open(avg, 'rb'))

    # Define color and name mappings for species
    color_dict = {
        'M_rdf': 'deeppink',
        'E1_rdf': 'yellowgreen',
        'E2_rdf': 'darkgreen',
        'P1_rdf': '#8C510A',
        'P2_rdf': '#5E918B',
        'P3_rdf': '#805D94'
    }
    name_dict = {
        'M_rdf': 'M',
        'E1_rdf': "5'",
        'E2_rdf': "3'",
        'P1_rdf': 'NONO',
        'P2_rdf': 'FUS',
        'P3_rdf': 'TDP43'
    }

    # Create subplots for each species
    if oldmodel:
        figure, axes = plt.subplots(1, len(color_dict) - 1, figsize=(24, 4))
    else:
        figure, axes = plt.subplots(1,len(color_dict), figsize=(24, 4))
    for seed in data.keys():
        i = 0
        for key in data[seed].keys():
            if key.endswith('rdf'):
                if oldmodel and key == 'P3_rdf':
                    continue  # Skip P3 for old models if specified

                # Plot average data
                axes[i].plot(
                    avgdata[key.split('_')[0] + '_bin'][-1],
                    avgdata[key][-1],
                    color='k',
                    linestyle='--',
                    linewidth=2,
                    label='Avg'
                )

                # Plot standard deviation if enabled
                if std:
                    # Calculate standard deviation
                    avg_values = np.array([data[s][key][-1] for s in data.keys()])
                    std_dev = np.std(avg_values, axis=0)

                    # Plot shaded region for standard deviation
                    axes[i].fill_between(
                        avgdata[key.split('_')[0] + '_bin'][-1],
                        np.maximum(avgdata[key][-1] - std_dev,0),
                        avgdata[key][-1] + std_dev,
                        color=color_dict[key],
                        alpha=0.1,
                        label='Std. Dev.'
                    )

                # Set plot title and limits
                axes[i].set_title(f'{name_dict[key]} Profile')
                axes[i].set_xlim(0, 1)

                # Ensure unique legend entries
                handles, labels = axes[i].get_legend_handles_labels()
                by_label = dict(zip(labels, handles))
                axes[i].legend(by_label.values(), by_label.keys(),loc='best')
                i += 1

    # Adjust layout and add labels
    figure.tight_layout(rect=[0.03, 0.05, 0.90, 0.90])
    figure.supxlabel('Normalized Distance from Center')
    figure.supylabel('Number Density')

    # Show the plot
    plt.show()

