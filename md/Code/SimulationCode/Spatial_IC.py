import numpy as np
import random


# Function to calculate Euclidean distance
def euclidean_distance(point1, point2):
    """
    Calculate the Euclidean distance between two points in 3D space.

    Parameters:
    point1 (tuple): The (x, y, z) coordinates of the first point.
    point2 (tuple): The (x, y, z) coordinates of the second point.

    Returns:
    float: The Euclidean distance between the two points.
    """
    return np.sqrt(np.sum((np.array(point1) - np.array(point2)) ** 2))

def check_overlaps(polymer_dict):
    """
    Check for overlaps in the polymer dictionary.

    Parameters:
    polymer_dict (dict): A dictionary where keys are polymer IDs and values are tuples of (x, y, z) coordinates of monomers.

    Returns:
    None
    """

    # List to store overlapping points
    overlapping_points = []

    # Check for overlaps between points in the same polymer
    for polymer_id, (x, y, z) in polymer_dict.items():
        points = list(zip(x, y, z))
        num_points = len(points)
        for i in range(num_points):
            for j in range(i + 1, num_points):
                if euclidean_distance(points[i], points[j]) < 1:
                    overlapping_points.append((polymer_id, points[i], points[j]))

    #If there are overlapping points, print the details
    if overlapping_points:
        print("Overlapping points or points with distance < 1 found:")
        for polymer_id, point1, point2 in overlapping_points:
            print(f"Polymer {polymer_id}: {point1} and {point2}")
    else:
        print("No overlapping points or points with distance < 1 found.")

def pre_calcs(NP,L,vcalc,walk_type,Nproteins=None,spacing=None,Lbox=None,vrna=None,seed=None,mdseed=None):
    """
    Pre-calculation function for polymer initialization in simulations. Seed and mdseed are logged for consistency in simulations.

    Parameters:
    - NP (int): Number of polymers
    - L (int): Length of each polymer (ie. number of monomers)
    - vcalc (str): Volume calculation method ('flexbox':constant volume concentration or 'fixedbox':contant box size)
    - walk_type (str): Type of walk ('multiplewalks' or 'line')
    - Nproteins (int, optional): Number of proteins
    - spacing (int, optional): Spacing between polymers
    - Lbox (float, optional): Length of the simulation box
    - vrna (float, optional): Volume fraction of RNA
    - seed (int, optional): Random seed for polymer generation
    - mdseed (int, optional): Random seed for molecular dynamics

    Returns:
    - Depending on vcalc and walk_type, returns different sets of parameters:
      - For 'flexbox': (size, seed, mdseed)
      - For 'fixedbox': (vprotein, vrna, seed, mdseed)
    """

    #generate seed if not provided
    if seed==None:
        seed=random.randint(0,2**32-1)
    if mdseed==None:
        mdseed=random.randint(0,65535)
    
    #generate multiple walks configuration, if it fails, generate a new seed and try again until a configuration with no overlaps made
    testdict=multiple_walks(NP,L,Lbox,spacing,seed)
    if testdict==None:
        while testdict==None:
            seed=random.randint(0,2**32-1)
            testdict=multiple_walks(NP,L,Lbox,spacing,seed)

    #Calculate properties for configuration based on walktype
    #Multiple walks: Different copies of NEAT1 in different configurations
    if walk_type=="multiplewalks":
        #if concentration is constant, calculate box size based on volume fraction
        if vcalc=='flexbox':
            Vsphere=4/3*np.pi*(0.5)**3
            Vrna=L*Vsphere
            Vbox=(NP*Vrna)/vrna
            size=Vbox**(1/3)
            return size,seed,mdseed
        
        #if box size is fixed, calculate volume fraction
        if vcalc=='fixedbox':
            Vprotein=(4/3)*np.pi*(0.25)**3 # diameter of protein is 0.5
            vprotein=(Nproteins*Vprotein)/(Lbox**3)
            Vsphere=4/3*np.pi*(0.5)**3
            Vrna=L*Vsphere
            vrna=(NP*Vrna)/(Lbox**3)
            return vprotein,vrna,seed,mdseed
    
    #Error check
    if Nproteins==None and Lbox==None and vrna==None:
        print("Error: Nproteins, Lbox and vrna need to be specified")
        return

def surrounding_check(x_p,y_p,z_p, b,polymer_check):
    """
    Check if the position (x_p, y_p, z_p) is completely surrounded by other occupied positions.
    This function checks the six possible moves (positive and negative directions along x, y, and z axes)
    from the given position (x_p, y_p, z_p) with step size `b`. If any of these positions are not occupied
    (i.e., not present in `polymer_check`), the function returns False, indicating that the position is not
    completely surrounded. If all possible moves are occupied, the function returns True.
    Parameters:
    x_p (int): The x-coordinate of the position to check.
    y_p (int): The y-coordinate of the position to check.
    z_p (int): The z-coordinate of the position to check.
    b (int): The step size for checking surrounding positions.
    polymer_check (set): A set of tuples representing occupied positions.
    Returns:
    bool: True if the position is completely surrounded by other occupied positions, False otherwise.
    """

    # Check the six possible moves from the current position
    possible_moves = [(b, 0, 0), (-b, 0, 0), (0, b, 0), (0, -b, 0), (0, 0, b), (0, 0, -b)]
    for move in possible_moves:
        new_x = x_p + move[0]
        new_y = y_p + move[1]
        new_z = z_p + move[2]

        # Check if the new position is not occupied 
        if (new_x, new_y, new_z) not in polymer_check: 
            # Found a move that is not occupied, so the position is not completely surrounded
            return False
    # All possible moves are occupied, so the position is completely surrounded
    return True

def multiple_walks(NP, L, box_length,spacing,seed):  
    """
    Generate multiple random walks for polymers within a simulation box.

    Parameters:
    - NP (int): Number of polymers
    - L (int): Length of each polymer (number of monomers)
    - box_length (float): Length of the simulation box
    - spacing (float): Region of replication. 
    - seed (int): Random seed for polymer generation. These are the seeds generated in the pre-calculation function.

    Returns:
    - dict: A dictionary where keys are polymer IDs and values are tuples of (x, y, z) coordinates of monomers
    """

    # Set random seed
    random.seed(seed)
    np.random.seed(seed)

    # Initialize dictionaries and variables
    polymer_dict={}
    b=1         
    polymer_check=set()  

    # Generate random walks for each polymer
    for polymer_id in range(NP):
        # Initialize polymer coordinates
        start_box_size=int(box_length/spacing) 
        x=np.zeros(L)
        y=np.zeros(L)
        z=np.zeros(L)
        fail_starts=set()
        start=False

        # Generate a random starting position for the polymer
        while not start:
            available_starts=[(x,y,z) for x in range(-start_box_size,start_box_size + 1) 
                             for y in range(-start_box_size,start_box_size +1)
                             for z in range(-start_box_size,start_box_size+1)
                             if (x,y,z) not in fail_starts]
            
            
            start_position=available_starts[np.random.randint(len(available_starts))]
            x[0],y[0],z[0]=start_position

            # Check if the starting position is already occupied
            if (x[0],y[0],z[0]) not in polymer_check:
                polymer_check.add((x[0],y[0],z[0]))
                start=True
            else:
                fail_starts.add(start_position)
        
        #If the starting position is valid, generate the rest of the polymer
        possible_steps=[(b,0,0),(-b,0,0),(0,b,0),(0,-b,0),(0,0,b),(0,0,-b)] 

        #Loop through the rest of the polymer
        for i in range(1, L):
            failed_steps=set()
            step_taken = False

            # Try to take a step until a valid step is taken
            while not step_taken:
                new_x, new_y, new_z = x[i - 1], y[i - 1], z[i - 1]
                available_steps=[step for step in possible_steps if step not in failed_steps] #this index is changing as fail steps occur
                if len(available_steps)==0:
                    return None
                # Randomly generate the step directions
                step_index=np.random.randint(len(available_steps))
                move=available_steps[step_index]
                # Calculate potential new position
                new_x += move[0]
                new_y += move[1]
                new_z += move[2]
                # Check if the new position is already occupied
                if ((new_x, new_y, new_z) not in polymer_check and
                    not surrounding_check(new_x, new_y, new_z, b, polymer_check) and
                    np.abs(new_x) < box_length / 2 and
                    np.abs(new_y) < box_length / 2 and
                    np.abs(new_z) < box_length / 2):
                    x[i], y[i], z[i] = new_x, new_y, new_z
                    polymer_check.add((new_x, new_y, new_z))
                    step_taken=True 
                else: 
                    failed_steps.add(move)
        #If the polymer is completed, add it to the dictionary
        polymer_dict[polymer_id]=(x,y,z)
    return polymer_dict

def initialize_polymers(NP, L,spacing,walk_type,seed,Lbox):
    """Initialize the polymer dictionary with random walk polymers"""

    if walk_type == "multiplewalks":
        return multiple_walks(NP, L, Lbox,spacing,seed)
    else:
        raise ValueError("Invalid walk type")