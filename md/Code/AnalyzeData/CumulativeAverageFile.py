import hoomd
import gsd.hoomd
import gsd
import freud
import numpy as np

class CumulativeAverage:
    """
    This class calculates the cumulative average of numbers, lists, numpy arrays, and dictionaries.
    
    Attributes:
        numbers (list): A list to store the numbers or arrays added.
        dict (dict): A dictionary to store the values when a dictionary is added.
        current_average (float, np.ndarray, dict): The current cumulative average.
        cumulative_averages (list): A list to store the cumulative averages after each addition.
    
        
    Methods:
        add_value: Add a new value to the cumulative average.
        get_cumulative_averages: Get the cumulative averages after each addition.
    """
    def __init__(self):
        self.numbers = []
        self.dict={}
        self.current_average = None
        self.cumulative_averages = []

    def add_value(self, new_value):
        if isinstance(new_value, (int,float)):
            self.numbers.append(new_value)
            sum_list=0
            if self.current_average is None:
                self.current_average = new_value
            else:
                for i in range(len(self.numbers)):
                    sum_list+=self.numbers[i]
                self.current_average = sum_list/(len(self.numbers))
                
        elif isinstance(new_value, (list,np.ndarray)):
            new_value=np.array(new_value)
            self.numbers.append(new_value)
            if self.current_average is None:
                self.current_average = new_value
            else:
                sum_list=0
                for i in range(len(self.numbers)):
                    sum_list+=self.numbers[i]
                self.current_average = sum_list/(len(self.numbers))

        elif isinstance(new_value, dict):
            for key, value in new_value.items():
                if key not in self.numbers:
                    self.dict[key] = []
                self.dict[key].append(value)
                if self.current_average is None:
                    self.current_average = {key: value}
                else:
                    sum_list=0
                    for i in range(len(self.dict[key])):
                        sum_list+=self.dict[key][i]
                    self.current_average[key] = sum_list/(len(self.dict[key]))
        self.cumulative_averages.append(self.current_average)
        return self.current_average

    def get_cumulative_averages(self):
        return self.cumulative_averages


class DataAverager:
    """
    This class handles the averaging of data from multiple files or multiple clusters for Maximum Intensity Plots.
    
    Attributes:
        data (dict): A dictionary to store the accumulated data.
        file_count (dict): A dictionary to store the count of files for each key.
    
    Methods:
        add_data: Add data to the dictionary.
        compute_average: Compute the average of the data.
    """
    def __init__(self):
        self.data = {}
        self.file_count = {}

    def add_data(self, data):

        if isinstance(data, dict):
            for key,d in data.items():
                if key not in self.data:
                    self.data[key] = np.array(d)
                    self.file_count[key] = 1
                else:
                    self.data[key] += np.array(d)
                    self.file_count[key] += 1
        elif isinstance(data, (list,np.ndarray)):
            if 'total' not in self.data:
                self.data['total'] = np.array(data)
                self.file_count['total'] = 1
            else:
                self.data['total'] += np.array(data)
                self.file_count['total'] += 1

    def compute_average(self):
        averaged_data = {key: value / self.file_count[key] for key, value in self.data.items()}
        return averaged_data


def average_rdf(dictionary):
    """
    This function calculates the average across multiple clusters for data. 

    Args:
        dictionary (dict): A dictionary where keys are cluster identifiers and values are dictionaries of RDF data.

    Returns:
        dict: A dictionary containing the averaged RDF data.
    """
    rdf_sum = {}
    rdf_sum = {}
    rdf_count = {}
    i = 0
    for cluster, data in dictionary.items():
        for p, rdf_values in data.items():
            if p not in rdf_sum:
                rdf_sum[p] = np.array(rdf_values)
                rdf_count[p] = 1
            else:
                rdf_sum[p] += np.array(rdf_values)
                rdf_count[p] += 1
        i += 1

    rdf_avg = {p: rdf_sum[p] / rdf_count[p] for p in rdf_sum}
    return rdf_avg


### Functions ###

def num_cluster_per_file(file):  
    """
    This function calculates the number of clusters per file, the average cluster mass, and the average radius of gyration.

    Args:
        file (str): The path to the GSD file.

    Returns:
        tuple: A tuple containing the cumulative averages of the number of clusters, cluster mass, cluster radius of gyration, and the simulation times.
    """
    RNA_TYPE_ID = {0, 1, 2}
    PROTEIN_TYPE_ID = {3, 4, 5}
    with gsd.hoomd.open(file, 'r') as f:
        num_clusters_avg = CumulativeAverage()
        cluster_mass_avg = CumulativeAverage()
        cluster_rg_avg = CumulativeAverage()
        times = []
        
        for i, frame in enumerate(f):
            times.append(frame.log['Simulation/timestep'][0])
            box = frame.configuration.box
            points = frame.particles.position
            types = frame.particles.typeid
            system = freud.AABBQuery(box, points)
            cl = freud.cluster.Cluster()
            cl.compute(system, neighbors={'r_max': 1.5})
            
            # Compute properties for valid clusters
            clp = freud.cluster.ClusterProperties()
            clp.compute(system, cl.cluster_idx)

            # Filter clusters to include only those with both RNA and protein
            valid_clusters = []
            for cluster_id in range(np.max(cl.cluster_idx) + 1):
                cluster_types = types[cl.cluster_idx == cluster_id]
                if any(t in RNA_TYPE_ID for t in cluster_types) and any(t in PROTEIN_TYPE_ID for t in cluster_types):
                    valid_clusters.append(cluster_id)
            num_clusters_avg.add_value(len(valid_clusters))

            # Compute properties for valid clusters
            valid_cluster_masses = np.mean([clp.cluster_masses[cluster_id] for cluster_id in valid_clusters])
            valid_cluster_rgs = np.mean([clp.radii_of_gyration[cluster_id] for cluster_id in valid_clusters])
            
            cluster_mass_avg.add_value(float(valid_cluster_masses))
            cluster_rg_avg.add_value(float(valid_cluster_rgs))
    
    return num_clusters_avg.get_cumulative_averages(), cluster_mass_avg.get_cumulative_averages(), cluster_rg_avg.get_cumulative_averages(), times


def normalize_bins(bin_centers, characteristic_length):
    """
    This function normalizes the bin centers by dividing them by the characteristic length.

    Args:
        bin_centers (np.ndarray): The bin centers to be normalized.
        characteristic_length (float): The characteristic length to normalize by.

    Returns:
        np.ndarray: The normalized bin centers.
    """
    return bin_centers / characteristic_length

def calc_PE(file):
    """
    This function calculates the potential energy (PE) from a GSD file.

    Args:
        file (str): The path to the GSD file.

    Returns:
        tuple: A tuple containing the simulation timesteps and the cumulative averages of the potential energy.
    """
    with gsd.hoomd.open(file, 'r') as f:
        if len(f)==0:
            return None
        PE=CumulativeAverage()
        timestep=[]
        for i,frame in enumerate(f):
            PE.add_value(frame.log['md/compute/ThermodynamicQuantities/potential_energy'][0])
            timestep.append(frame.log['Simulation/timestep'][0])
        # Calculate running average
        PE_avg = PE.get_cumulative_averages()
        timestep=timestep
    return timestep,PE_avg
    
def MIP(file,valid_clusters,num_frames,axis='x'):
    """
    This function calculates the maximum intensity projection (MIP) graph for a given cluster along a specified axis.

    Args:
        file (str): The path to the GSD file.
        valid_clusters (list): A list of valid cluster IDs.
        num_frames (int): The number of frames to consider.
        axis (str): The axis along which to calculate the MIP graph ('x', 'y', or 'z').

    Returns:
        dict: A dictionary containing the averaged MIP graphs for each species and the total density.
    """

    #Define particle identities and size for each and together
    dictionary = {'M': 0, 'E1': 1, 'E2': 2, 'P1': 3, 'P2': 4, 'P3': 5}
    volportein=(4/3)*np.pi*(0.25**3)
    volrna=(4/3)*np.pi*(0.5**3)
    total_volume=np.mean(volportein+volrna)
    #Calculate MIP for each cluster
    cluster_mip=DataAverager()
    matricies={}
    for cluster_id in valid_clusters:
        #calculate the MIP for each frame of cluster_id
        clusterframe=CumulativeAverage()
        with gsd.hoomd.open(file) as f:
            for frame in f[-num_frames:]:
                box = frame.configuration.box
                points = frame.particles.position
                typeids = frame.particles.typeid
                
                # Compute local density for the entire system
                ld = freud.density.LocalDensity(r_max=1, diameter=1.0)
                ld.compute((box, points))
                
                # Compute clusters
                cl = freud.cluster.Cluster()
                cl.compute((box, points), neighbors={'r_max': 1.5})
                
                # Compute cluster properties
                clp = freud.cluster.ClusterProperties()
                clp.compute((box, points), cl.cluster_idx)
                
                # Get the center and radius of gyration of the specified cluster
                centers = clp.centers[cluster_id]
                rg = clp.radii_of_gyration[cluster_id]
                
                # Define the grid size based on the radius of gyration
                grid_size_factor = 2  # Factor to make the grid slightly wider than the cluster
                grid_size = int(rg * grid_size_factor)
                
                # Define the grid
                x_min, x_max = centers[0] - grid_size, centers[0] + grid_size
                y_min, y_max = centers[1] - grid_size, centers[1] + grid_size
                z_min, z_max = centers[2] - grid_size, centers[2] + grid_size
                grid_resolution = 30  # Define the grid resolution
                x_bins = np.linspace(x_min, x_max, grid_resolution)
                y_bins = np.linspace(y_min, y_max, grid_resolution)
                z_bins = np.linspace(z_min, z_max, grid_resolution)
                max_density_grid = np.zeros((grid_resolution - 1, grid_resolution - 1))
                num_slices = 300


                #Define the axis for MIP
                if axis == 'x':
                    slices = np.linspace(x_min, x_max, num_slices + 1)
                    index=int(0)
                    bin1=y_bins
                    bin2=z_bins
                    in1=int(1)
                    in2=int(2)
                elif axis == 'y':
                    slices = np.linspace(y_min, y_max, num_slices + 1)
                    index=int(1)
                    bin1=x_bins
                    bin2=z_bins
                    in1=int(0)
                    in2=int(2)
                elif axis == 'z':
                    slices = np.linspace(z_min, z_max, num_slices + 1)
                    index=int(2)
                    bin1=x_bins
                    bin2=y_bins
                    in1,in2=int(0),int(1)

                for i in range(num_slices):
                    start = slices[i]
                    end = slices[i + 1]
                    slice_atoms = points[(points[:, index] >= start) & (points[:, index] < end)]
                    slice_density = ld.density[(points[:, index] >= start) & (points[:, index] < end)]
                    if len(slice_atoms) == 0:
                        continue
                    heatmap, _, _ = np.histogram2d(slice_atoms[:, in1], slice_atoms[:, in2], bins=[bin1, bin2], weights=slice_density*total_volume)
                    max_density_grid = np.maximum(max_density_grid, heatmap)
                    max_density_grid[max_density_grid < 0] = np.nan
                matricies['total']=max_density_grid

                # Calculate the MIP for each species through the slices of simulation space
                for species in dictionary:
                    max_density_grid = np.zeros((grid_resolution - 1, grid_resolution - 1))

                    #Define local density for each species
                    species_pos = points[typeids == dictionary[species]]
                    species_density = ld.density[typeids == dictionary[species]]

                    #Define volume for each species
                    if species=='M' or species=='E1' or species=='E2':
                        volume=volrna
                    elif species=='P1' or species=='P2' or species=='P3':
                        volume=volportein

                    #Rotate through slices to calculate MIP
                    for i in range(num_slices):
                        start = slices[i]
                        end = slices[i + 1]
                        slice_atoms = species_pos[(species_pos[:, index] >= start) & (species_pos[:, index] < end)]
                        slice_density = species_density[(species_pos[:, index] >= start) & (species_pos[:, index] < end)]
                        if len(slice_atoms) == 0:
                            continue
                        heatmap, _, _ = np.histogram2d(slice_atoms[:, in1], slice_atoms[:, in2], bins=[bin1, bin2], weights=slice_density*volume)
                        max_density_grid = np.maximum(max_density_grid, heatmap)
                        max_density_grid[max_density_grid < 0] = np.nan
                    
                    #store the MIP for each species
                    matricies[str(species)]=max_density_grid
                #store the MIP for cluster_id at frame f
                clusterframe.add_value(matricies)
        #store the MIP for cluster_id across all frames
        cluster_mip.add_data(clusterframe.get_cumulative_averages()[0])

    #average the MIP for all clusters
    return cluster_mip.compute_average()


def calc_density(file,valid_clusters,num_frames):
    """
    This function calculates the intensity of different species within clusters over multiple frames.

    Args:
        file (str): The path to the GSD file.
        valid_clusters (list): A list of valid cluster IDs.
        num_frames (int): The number of frames to consider.

    Returns:
        dict: A dictionary containing the averaged radial density and bin centers for each species.
    """
    clusters_den={}
    # Loop through each cluster_id
    for cluster_id in valid_clusters:   
        with gsd.hoomd.open(file, 'r') as f:
            M_avg=CumulativeAverage()
            M_bins=CumulativeAverage()
            E1_avg=CumulativeAverage()
            E1_bins=CumulativeAverage()
            E2_avg=CumulativeAverage()
            E2_bins=CumulativeAverage()
            P1_avg=CumulativeAverage()
            P1_bins=CumulativeAverage()
            P2_avg=CumulativeAverage()
            P2_bins=CumulativeAverage()
            P3_avg=CumulativeAverage()
            P3_bins=CumulativeAverage()

            # Loop through each frame
            for frame in f[-num_frames:]:
                box = frame.configuration.box
                points = frame.particles.position
                particle_types = frame.particles.typeid
                

                #find clusters
                total_system= freud.AABBQuery(box, points)
                cl = freud.cluster.Cluster()
                cl.compute(total_system, neighbors={'r_max': 1.5})
                clp = freud.cluster.ClusterProperties()
                clp.compute(total_system, cl.cluster_idx)
                center=clp.centers[cluster_id]


                #define cluster rg and characterisitic length
                rg=clp.radii_of_gyration[cluster_id]
                characteristic_length = int(rg*2)

                #Specific clusters in frame particles positions and types
                cluster_indices = np.where(cl.cluster_idx == cluster_id)[0]
                cluster_particles = points[cluster_indices]
                num_bins = 10
                bin_edges = np.linspace(0, characteristic_length, num_bins + 1)
                bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
                

                # Filter cluster particles by type
                M_cluster = cluster_particles[particle_types[cluster_indices] == 0]
                E1_cluster = cluster_particles[particle_types[cluster_indices] == 1]
                E2_cluster = cluster_particles[particle_types[cluster_indices] == 2]
                P1_cluster = cluster_particles[particle_types[cluster_indices] == 3]
                P2_cluster = cluster_particles[particle_types[cluster_indices] == 4]
                P3_cluster = cluster_particles[particle_types[cluster_indices] == 5]
                species=[M_cluster,E1_cluster,E2_cluster,P1_cluster,P2_cluster,P3_cluster]
                bins=[M_bins,E1_bins,E2_bins,P1_bins,P2_bins,P3_bins]
                rdfs=[M_avg,E1_avg,E2_avg,P1_avg,P2_avg,P3_avg]
                densities = {}
                
                # Calculate the radial distribution function (RDF) for each species
                for sp,bin,avg in zip(species,bins,rdfs):
                    distances = np.linalg.norm(sp - center, axis=1)  
                    # Bin the distances
                    hist, _ = np.histogram(distances, bins=bin_edges)

                    # Normalize the histogram
                    shell_volumes = (4/3) * np.pi * (bin_edges[1:] ** 3 - bin_edges[:-1] ** 3)
                    rdf = hist/shell_volumes



                    avg.add_value(rdf)
                    bin.add_value(normalize_bins(bin_centers, characteristic_length))
                densities['M_rdf'] = M_avg.get_cumulative_averages()
                densities['M_bin'] = M_bins.get_cumulative_averages()
                densities['E1_rdf'] = E1_avg.get_cumulative_averages()
                densities['E1_bin'] = E1_bins.get_cumulative_averages()
                densities['E2_rdf'] = E2_avg.get_cumulative_averages()
                densities['E2_bin'] = E2_bins.get_cumulative_averages()
                densities['P1_rdf'] = P1_avg.get_cumulative_averages()
                densities['P1_bin'] = P1_bins.get_cumulative_averages()
                densities['P2_rdf'] = P2_avg.get_cumulative_averages()
                densities['P2_bin'] = P2_bins.get_cumulative_averages()
                densities['P3_rdf'] = P3_avg.get_cumulative_averages()
                densities['P3_bin'] = P3_bins.get_cumulative_averages()
            
            clusters_den[cluster_id]=densities

    cluster_int_avg=average_rdf(clusters_den)


    density_rdf={
        'M_rdf':cluster_int_avg['M_rdf'],
        'M_bin':cluster_int_avg['M_bin'],
        'E1_rdf':cluster_int_avg['E1_rdf'],
        'E1_bin':cluster_int_avg['E1_bin'],
        'E2_rdf':cluster_int_avg['E2_rdf'],
        'E2_bin':cluster_int_avg['E2_bin'],
        'P1_rdf':cluster_int_avg['P1_rdf'],
        'P1_bin':cluster_int_avg['P1_bin'],
        'P2_rdf':cluster_int_avg['P2_rdf'],
        'P2_bin':cluster_int_avg['P2_bin'],
        'P3_rdf':cluster_int_avg['P3_rdf'],
        'P3_bin':cluster_int_avg['P3_bin'],
    }
    

    return density_rdf


def calc_density_along_rna(file,valid_clusters,num_frames):
    """
    This function calculates the density of RNA segments (first, middle, last) within clusters over multiple frames.

    Args:
        file (str): The path to the GSD file.
        valid_clusters (list): A list of valid cluster IDs.
        num_frames (int): The number of frames to consider.

    Returns:
        dict: A dictionary containing the averaged RDFs and bin centers for each RNA segment (first, middle, last).
    """
    #Loop through valid clusters
    cluster_rna_int={}
    for cluster_id in valid_clusters:
        with gsd.hoomd.open(file, 'r') as f:
            first_avg=CumulativeAverage()
            first_bins=CumulativeAverage()
            middle_avg=CumulativeAverage()
            middle_bins=CumulativeAverage()
            last_avg=CumulativeAverage()
            last_bins=CumulativeAverage()

            # Loop through each frame
            i=0
            for frame in f[-num_frames:]:

                #frame box
                box = frame.configuration.box
                points = frame.particles.position
                particle_types = frame.particles.typeid
                total_system= freud.AABBQuery(box, points)
                
                #find clusters
                cl = freud.cluster.Cluster()
                cl.compute(total_system, neighbors={'r_max': 1.5})
                clp = freud.cluster.ClusterProperties()
                clp.compute(total_system, cl.cluster_idx)
                center=clp.centers[cluster_id]
                rg=clp.radii_of_gyration[cluster_id]

                # Calculate characteristic length scale for binning
                characteristic_length = int(2*rg) 
                num_bins = 10
                bin_edges = np.linspace(0, characteristic_length, num_bins + 1)
                bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

                #find cluster particles for specific cluster_id
                cluster_indices = np.where(cl.cluster_idx == cluster_id)[0]
                cluster_particles = points[cluster_indices]
                cluster_types = particle_types[cluster_indices]

                # Filter cluster particles by RNA species
                rna_species = [0, 1, 2]
                rna_indices = np.isin(cluster_types, rna_species)
                rna_particles = cluster_particles[rna_indices]
                num_poly=int(len(rna_particles)//45)

                # Split RNA particles into first, middle, and last segments
                first15=[]
                middle15=[]
                last15=[]
                for k in range(num_poly):
                    polymer = rna_particles[k * 45:(k + 1) * 45]
                    first15.extend(polymer[:15])
                    middle15.extend(polymer[15:30])
                    last15.extend(polymer[30:45])
                species=[first15,middle15,last15]
                rdfs=[first_avg,middle_avg,last_avg]
                bins=[first_bins,middle_bins,last_bins]
                densities = {}
                for sp,avg,bin in zip(species,rdfs,bins):
                    #calculate distance from center
                    distances = np.linalg.norm(sp - center, axis=1)
                    # Bin the distances
                    hist, _ = np.histogram(distances, bins=bin_edges)
                    # Normalize the histogram
                    shell_volumes = (4/3) * np.pi * (bin_edges[1:] ** 3 - bin_edges[:-1] ** 3)
                    rdf = hist/shell_volumes 
                    avg.add_value(rdf)
                    bin.add_value(normalize_bins(bin_centers, characteristic_length))
                #store the MIP for cluster_id at frame f
                densities['first_rdf'] = first_avg.get_cumulative_averages()
                densities['first_bin'] = first_bins.get_cumulative_averages()
                densities['middle_rdf'] = middle_avg.get_cumulative_averages()
                densities['middle_bin'] = middle_bins.get_cumulative_averages()
                densities['last_rdf'] = last_avg.get_cumulative_averages()
                densities['last_bin'] = last_bins.get_cumulative_averages()
                i+=1
            #store across clusters
            cluster_rna_int[cluster_id]=densities
    #average across all clusters
    cluster_rna_int_avg=average_rdf(cluster_rna_int)

    density_rna_rdf={
        'first_rdf':cluster_rna_int_avg['first_rdf'],
        'first_bin':cluster_rna_int_avg['first_bin'],
        'middle_rdf':cluster_rna_int_avg['middle_rdf'],
        'middle_bin':cluster_rna_int_avg['middle_bin'],
        'last_rdf':cluster_rna_int_avg['last_rdf'],
        'last_bin':cluster_rna_int_avg['last_bin'],
    }
        
    return density_rna_rdf




