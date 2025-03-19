import hoomd
import gsd
import gsd.hoomd
import os
import random
import numpy as np
import Spatial_IC as ic
import DataMangement as dm2p_ri

def create_sys(inputpath):
    """Create the initial system configuration with NEAT1 RNA and proteins. 
    Inputs:
    - inputpath: path to the input file from JobInputs folder
    Outputs:
    - .gsd file with the initial configuration of the system
    """

    print("Creating system")

    #Read input parameters
    params=dm2p_ri.input_params(inputpath)
    NP=int(params['Npolymer'])
    L=int(params['L'])
    Nproteins=int(params['Nproteins'])
    LE=int(params['LE'])
    LM=int(params['LM'])
    spacing=params['spacing']
    seed=int(params['seed'])
    vP1=int(params['vP1'])
    vP2=int(params['vP2'])
    vP3=int(params['vP3'])
    walk_type=params['WalkType']
    Lbox=int(params['Lbox'])


    #making polymers in initialized configuration
    polymer_dict=ic.initialize_polymers(NP,L,spacing,walk_type,seed,Lbox)

    #set seed for random
    np.random.seed(seed)
    random.seed(seed)
    
    ###Create frame and initialize particles
    frame = gsd.hoomd.Frame()
    frame.particles.N = L*NP + Nproteins
    frame.particles.mass = np.ones(frame.particles.N)
    frame.particles.diameter = np.ones(frame.particles.N)
    frame.configuration.box=[Lbox,Lbox,Lbox,0,0,0]

    # Position RNA monomers in the box
    positions=[]
    for polymer_id in polymer_dict:
        x,y,z=polymer_dict[polymer_id]
        positions.extend([(xi,yi,zi) for xi,yi,zi in zip(x,y,z)])
    positions=np.array(positions)
    frame.particles.position=positions

    # Function to check for overlaps
    def check_overlap(new_pos, existing_positions, min_distance):
        for pos in existing_positions:
            if np.linalg.norm(new_pos - pos) < min_distance:
                return True
        return False

    # Randomly add proteins without interfering with RNA
    protein_positions = []
    min_distance_RNA = 0.75  # Minimum distance to avoid overlap
    min_distance_protein = 0.5  # Minimum distance between proteins

    # Add proteins to the box
    while len(protein_positions) < Nproteins:
        new_pos = np.random.uniform(low=-Lbox/2, high=Lbox/2, size=3)
        if not check_overlap(new_pos, positions, min_distance_RNA) and not check_overlap(new_pos, protein_positions, min_distance_protein):
            protein_positions.append(new_pos)
    protein_positions = np.array(protein_positions)
    frame.particles.position = np.vstack((positions, protein_positions))

    # Update diameters and masses for proteins
    frame.particles.diameter[-Nproteins:] = 0.5
    
    ### Define Particle Types
    frame.particles.types = ['A', 'B', 'C','P1','P2','P3']
    frame.particles.typeid = [0] * frame.particles.N #initialize all particles as A

    # Assign particle types for proteins and shuffle
    proteins_types=[3]*vP1+[4]*vP2+[5]*vP3
    random.shuffle(proteins_types)

    # Assign particle types for RNA
    for i in range(NP):
        typeids = [1] * LE + [0] * LM + [2] * LE
        frame.particles.typeid[i * L:(i + 1) * L] = typeids
    
    # Assign particle types for proteins
    for k, protein_type in zip(range(-Nproteins, 0), proteins_types):
        frame.particles.typeid[k] = protein_type
    

    ### Define Bonds and Forces
    #define harmonic bonds (connected particles)
    frame.bonds.N = (L-1)*NP
    frame.bonds.types=['A-A']
    frame.bonds.typeid=[0]*(L-1)*NP

    # Define bonds and forces
    frame.bonds.N = (L - 1) * NP
    frame.bonds.types = ['A-A']
    frame.bonds.typeid = [0] * frame.bonds.N
    bond_groups = []
    for n in range(NP):
        start_index = n * L
        for i in range(start_index, start_index + L - 1):
            bond_groups.append([i, i + 1])
    frame.bonds.group = bond_groups
    
    #create initial state
    directory,base=os.path.split(params['InitialStateFilePath'])
    os.makedirs(directory,exist_ok=True)

    #save intitial state as a .gds
    with gsd.hoomd.open(name=params['InitialStateFilePath'], mode='w') as f:
        f.append(frame)
        
def run_traj(inputpath):
    """For an initial state .gsd file add harmonic and LJ forces. Run the simulation trajectory with two different logging periods.
    Inputs:
    - inputpath: path to the input file from JobInputs folder
    Outputs:
    - .gsd file with the trajectory of the system
    """
    print("Running trajectory")

    #Read input parameters
    params=dm2p_ri.input_params(inputpath)

    ##Define Forces
    #harmonic bonds  
    harmonic = hoomd.md.bond.Harmonic()
    harmonic.params['A-A']=dict(k=params['hk'],r0=params['r0'])

   #lennard jones
   #RNA-RNA
    lj = hoomd.md.pair.LJ(nlist=hoomd.md.nlist.Cell(buffer=0.4),default_r_cut=params['r_cut'],mode='shift')
    lj.params[('A','A')] = dict(epsilon=params['Mep'],sigma=params['Msigma'])
    lj.params[('B','B')] = dict(epsilon=params['Eep'],sigma=params['Esigma'])
    lj.params[('C','C')] = dict(epsilon=params['Eep'],sigma=params['Esigma'])
    lj.params[('A','B')] = dict(epsilon=params['EMep'],sigma=params['EMsigma'])
    lj.params[('A','C')] = dict(epsilon=params['EMep'],sigma=params['EMsigma'])
    lj.params[('B','C')] = dict(epsilon=params['Eep'],sigma=params['Esigma'])
    #RNA-P1
    lj.params[('A','P1')]=dict(epsilon=params['P1Mep'],sigma=params['PMsigma'])
    lj.params[('B','P1')]=dict(epsilon=params['P1E1ep'],sigma=params['PEsigma'])
    lj.params[('C','P1')]=dict(epsilon=params['P1E2ep'],sigma=params['PEsigma'])
    #RNA-P2
    lj.params[('A','P2')]=dict(epsilon=params['P2Mep'],sigma=params['PMsigma'])
    lj.params[('B','P2')]=dict(epsilon=params['P2E1ep'],sigma=params['PEsigma'])
    lj.params[('C','P2')]=dict(epsilon=params['P2E2ep'],sigma=params['PEsigma'])
    #RNA-P3
    lj.params[('A','P3')]=dict(epsilon=params['P3Mep'],sigma=params['PMsigma'])
    lj.params[('B','P3')]=dict(epsilon=params['P3E1ep'],sigma=params['PEsigma'])
    lj.params[('C','P3')]=dict(epsilon=params['P3E2ep'],sigma=params['PEsigma'])
    #protein-proteins
    lj.params[('P1','P1')]=dict(epsilon=params['P1ep'],sigma=params['Psigma'])
    lj.params[('P2','P2')]=dict(epsilon=params['P2ep'],sigma=params['Psigma'])
    lj.params[('P3','P3')]=dict(epsilon=params['P3ep'],sigma=params['Psigma'])
    lj.params[('P1','P2')]=dict(epsilon=params['P12ep'],sigma=params['Psigma'])
    lj.params[('P1','P3')]=dict(epsilon=params['P13ep'],sigma=params['Psigma'])
    lj.params[('P2','P3')]=dict(epsilon=params['P23ep'],sigma=params['Psigma'])

    force=[harmonic,lj]

    #set up sim
    device=hoomd.device.GPU()
    #set seed: if usemdseed is 1, use the mdseed in input file, otherwise generate a new random seed
    if int(params['usemdseed'])==1:
        mdseed=int(params['mdseed'])
    else:
        mdseed=random.randint(0,2**32-1)

    #create sim
    sim=hoomd.Simulation(device=device,seed=mdseed)
    sim.create_state_from_gsd(filename=params['InitialStateFilePath'])
    langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(),kT=1.0)
    integrator=hoomd.md.Integrator(dt=params['dt'],methods=[langevin],forces=force)
    sim.operations.integrator=integrator

    ##Logging##
    #add logger
    thermo=hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)
    logger=hoomd.logging.Logger()
    logger.add(sim,quantities=['timestep','tps','walltime'])
    logger.add(thermo)

    #make traj file
    outpath=params['TrajectoryFilePath']
    os.makedirs(os.path.dirname(outpath),exist_ok=True)
    
    #log gsd
    ##write initial frame
    gsd_writer_initial=hoomd.write.GSD(trigger=hoomd.trigger.Periodic(int(params['initialperiod'])),
            filename=outpath,
            filter=hoomd.filter.All(),
            mode='ab',
            logger=logger)

    sim.operations.writers.append(gsd_writer_initial)
    sim.run(params['firsthalf'])
    gsd_writer_initial.flush()
    
    ##write final frame
    sim.operations.writers.remove(gsd_writer_initial)
    gsd_writer_final=hoomd.write.GSD(trigger=hoomd.trigger.Periodic(int(params['finalperiod'])),
            filename=outpath,
            filter=hoomd.filter.All(),
            mode='ab',
            logger=logger)

    sim.operations.writers.append(gsd_writer_final)
    sim.run(params['secondhalf'])
    gsd_writer_final.flush()

    
def run_sim(inputpath):
    """Run the simulation based on the input file.
    Inputs:
    - inputpath: path to the input file from JobInputs folder
    Outputs:
    - .gsd file with the trajectory of the system
    """
    create_sys(inputpath)
    run_traj(inputpath)
    
    
    
    
