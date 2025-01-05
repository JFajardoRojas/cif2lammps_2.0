import networkx as nx
import numpy as np
import atomic_data
import operator
from itertools import groupby, combinations, permutations
from random import randint
import small_molecule_constants
from cif2system import PBC3DF_sym
import write_molecule_files as WMF
from ase import Atom, Atoms
from ase import neighborlist
from ase.geometry import get_distances
from ase.io import read

#Just to test the graphs
import matplotlib.pyplot as plt

mass_key = atomic_data.mass_key

def add_small_molecules(FF, ff_string):
    
    if ff_string == 'TraPPE':
        SM_constants = small_molecule_constants.TraPPE
        FF.pair_data['special_bonds'] = 'lj 0.0 0.0 1.0 coul 0.0 0.0 0.0'   #Including to meet the tests done on MOF+Solvent

    elif ff_string == 'TIP4P_2005_long':
        SM_constants = small_molecule_constants.TIP4P_2005_long
        FF.pair_data['special_bonds'] = 'lj 0.0 0.0 1.0 coul 0.0 0.0 0.0'

    elif ff_string == 'TIP4P_cutoff':
        SM_constants = small_molecule_constants.TIP4P_cutoff
        FF.pair_data['special_bonds'] = 'lj/coul 0.0 0.0 1.0'

    elif ff_string == 'TIP4P_2005_cutoff':
        SM_constants = small_molecule_constants.TIP4P_2005_cutoff
        FF.pair_data['special_bonds'] = 'lj/coul 0.0 0.0 1.0'

    elif ff_string == 'Ions':
        SM_constants = small_molecule_constants.Ions
        FF.pair_data['special_bonds'] = 'lj/coul 0.0 0.0 1.0'

    elif ff_string == 'UFF':
        SM_constants = small_molecule_constants.UFF
        FF.pair_data['special_bonds'] = 'lj 0.0 0.0 1.0 coul 0.0 0.0 0.0'


######## insert more force fields here if needed

    else:
        raise ValueError('the small molecule force field', ff_string, 'is not defined')

    SG = FF.system['graph']
    SMG = FF.system['SM_graph']

    # This is written based on the resulted CIF from the 2retrieve_DMF_loaded_structures.py    

    if len(SMG.nodes()) > 0 and len(SMG.edges()) == 0:   
        
        print('There are no small molecule bonds in the CIF, calculating based on covalent radii...')
        
        atoms = Atoms()

        offset = min(SMG.nodes())

        for node,data in SMG.nodes(data=True):
            atoms.append(Atom(data['element_symbol'], data['cartesian_position']))

        atoms.set_cell(FF.system['box'])
        unit_cell = atoms.get_cell()
        cutoffs = neighborlist.natural_cutoffs(atoms)
#        print(cutoffs)
        NL = neighborlist.NewPrimitiveNeighborList(cutoffs, use_scaled_positions=False, self_interaction=False, skin=0.10) # shorten the cutoff a bit
        NL.build([True, True, True], unit_cell, atoms.get_positions())

        for i in atoms:
            
            nbors = NL.get_neighbors(i.index)[0]

            for j in nbors:

                bond_length = get_distances(i.position, p2=atoms[j].position, cell=unit_cell, pbc=[True,True,True])
                bond_length = np.round(bond_length[1][0][0], 3)
                SMG.add_edge(i.index + offset, j + offset, bond_length=bond_length, bond_order='1.0', bond_type='S')

        NMOL = len(list(nx.connected_components(SMG)))
        print(NMOL, 'small molecules were recovered after bond calculation')

    mol_flag = 1
    max_ind = FF.system['max_ind']
    index = max_ind                


    box = FF.system['box']
    a,b,c,alpha,beta,gamma = box
    pi = np.pi
    ax = a
    ay = 0.0
    az = 0.0
    bx = b * np.cos(gamma * pi / 180.0)
    by = b * np.sin(gamma * pi / 180.0)
    bz = 0.0
    cx = c * np.cos(beta * pi / 180.0)
    cy = (c * b * np.cos(alpha * pi /180.0) - bx * cx) / by
    cz = (c ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
    unit_cell = np.asarray([[ax,ay,az],[bx,by,bz],[cx,cy,cz]]).T
    inv_unit_cell = np.linalg.inv(unit_cell)


    add_nodes = []
    add_edges = []
    comps = []

    for comp in nx.connected_components(SMG):

        mol_flag += 1
        comp = sorted(list(comp))
        ID_string = sorted([SMG.nodes[n]['element_symbol'] for n in comp])
        ID_string = [(key, len(list(group))) for key, group in groupby(ID_string)]
        ID_string = ''.join([str(e) for c in ID_string for e in c])
        comps.append(ID_string)

#        print(comps)

        for n in comp:

            data = SMG.nodes[n]
#            print(data)                # Checking the presence of charge in the data of the Molecule Graph

            SMG.nodes[n]['mol_flag'] = str(mol_flag)
#            print(SMG.nodes[n])
#            print(ID_string)

            if ID_string == 'H2O1':     #Recognition of the Water molecule by Graph Theory
                SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + '_w'

            elif ID_string == 'C1H1O1': #Recognition of the Methanol molecule by Graph Theory
                SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + '_CH3OH'

            elif ID_string == 'C6':  #Recognition of Hexane molecule by GT, Here Ch3 is Pd and Ch2 is Ga

                if SMG.nodes[n]['index'] == comp[0]:
                    SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + 'H3_HEX'   #This is written to recognize extreme elements of the graph as CH3s
                
                elif SMG.nodes[n]['index'] == comp[1]:
                    SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + 'H2_HEX'   #This is written to recognize extreme elements of the graph as CH2s

                elif SMG.nodes[n]['index'] == comp[2]:
                    SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + 'H2_HEX'   #This is written to recognize extreme elements of the graph as CH2s

                elif SMG.nodes[n]['index'] == comp[3]:
                    SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + 'H2_HEX'   #This is written to recognize extreme elements of the graph as CH2s

                elif SMG.nodes[n]['index'] == comp[4]:
                    SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + 'H2_HEX'   #This is written to recognize extreme elements of the graph as CH2s

                elif SMG.nodes[n]['index'] == comp[5]:
                    SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + 'H3_HEX'   #This is written to recognize extreme elements of the graph as CH3s

            
            elif ID_string == 'C3H7N1O1':  # This is the string representation of DMF
                # I am using the charge of each DMF atom to distinguish them. 
                if SMG.nodes[n]['charge'] == 0.0866:    # Charge for H in CH3
                    SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + '_CH3_dmf'
#                    print(data)
                elif SMG.nodes[n]['charge'] == 0.1234:    # Charge for H in COH
                    SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + '_COH_dmf'
#                    print(data)
                elif SMG.nodes[n]['charge'] == 0.1429:    # Charge for C in CH3
                    SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + '_CH3_dmf'
#                    print(data)
                elif SMG.nodes[n]['charge'] == 0.5924:    # Charge for C in COH
                    SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + '_COH_dmf'
#                    print(data)
                else:                                    # Nitrogen and Oxygen can be identified just using the name.
                    SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + '_dmf' 
#                    print(data)
            else:
                SMG.nodes[n]['force_field_type'] = SMG.nodes[n]['element_symbol'] + '_' + ID_string
#                print(data)

        # add COM sites where relevant, extend this list as new types are added
        if ID_string in ('O2', 'N2'):
            coords = []
            anchor = SMG.nodes[comp[0]]['fractional_position']

            for n in comp:

                data = SMG.nodes[n]
                data['mol_flag'] = str(mol_flag)
                fcoord = data['fractional_position']
                mic = PBC3DF_sym(fcoord, anchor)
                fcoord += mic[1]
                ccoord = np.dot(unit_cell, fcoord)
                coords.append(ccoord)

            ccom = np.average(coords, axis=0)
            fcom = np.dot(inv_unit_cell, ccom)
            index += 1

            if ID_string == 'O2':
                fft = 'O_com'
            elif ID_string == 'N2':
                fft = 'N_com'

            ndata =  {'element_symbol':'NA', 'mol_flag':mol_flag, 'index':index, 'force_field_type':fft, 'cartesian_position':ccom, 'fractional_position':fcom, 'charge':0.0, 'replication':np.array([0.0,0.0,0.0]), 'duplicated_version_of':None}
            edata =  {'sym_code':None, 'length':None, 'bond_type':None}

            add_nodes.append([index, ndata])
            add_edges.extend([(index, comp[0], edata), (index, comp[1], edata)])


    for n, data in add_nodes:
        SMG.add_node(n, **data)
    for e0, e1, data in add_edges:
        SMG.add_edge(e0, e1, **data)


    ntypes = max([FF.atom_types[ty] for ty in FF.atom_types])
#    print(ntypes)
    maxatomtype_wsm = max([FF.atom_types[ty] for ty in FF.atom_types])
#    print(maxatomtype_wsm)
    maxbondtype_wsm = max([bty for bty in FF.bond_data['params']])
    maxangletype_wsm = max([aty for aty in FF.angle_data['params']])

    nbonds = max([i for i in FF.bond_data['params']])
    nangles = max([i for i in FF.angle_data['params']])
  
    try:
        ndihedrals = max([i for i in FF.dihedral_data['params']])
#        print(ndihedrals)
    except ValueError:
        ndihedrals = 0

    try:
        nimpropers = max([i for i in FF.improper_data['params']])
#        print(nimpropers)
    except ValueError:
        nimpropers = 0

# Here is where the data for the detected molecules is included
    new_bond_types = {}
    new_angle_types = {}
    new_dihedral_types = {}
    new_improper_types = {}

    for subG, ID_string in zip([SMG.subgraph(c).copy() for c in nx.connected_components(SMG)], comps):

        constants = SM_constants[ID_string]
        #print(constants)


######### add new atom types

        for name,data in sorted(subG.nodes(data=True), key=lambda x:x[0]):

            fft = data['force_field_type']

            chg = constants['pair']['charges'][fft]
            data['charge'] = chg
            #print(fft, data['charge'])

            SG.add_node(name, **data)

            try: 
                FF.atom_types[fft] += 0
            except KeyError:

                ntypes += 1
                FF.atom_types[fft] = ntypes
                style = constants['pair']['style']
                vdW = constants['pair']['vdW'][fft]
#                print(vdW, fft)
                FF.pair_data['params'][FF.atom_types[fft]] = (style, vdW[0], vdW[1])
                FF.pair_data['comments'][FF.atom_types[fft]] = [fft, fft]
#                print(FF.atom_masses)
                
                if data['element_symbol'] == 'C' and ID_string == 'C1H1O1': #Modification made to recognize Methanol with pseudoatom CH3
                    print('##########')
                    print('Mass of CH3 for recognized molecule')
                    print('If you are not trying to include Methanol (CH3OH), you need to modify the mass assigment accordingly...')
                    FF.atom_masses[fft] = '15.03422'
                    print('##########')
                    #print(ID_string)

                elif ID_string == 'C6' and fft == 'CH3_HEX': #Modification made to recognize Ch3 in Hexane
                    FF.atom_masses[fft] = '15.03422'

                    print('##########')
                    print('Mass of CH3 for recognized molecule = ' +str(FF.atom_masses[fft]))
                    print('If you are not trying to include Hexane (C6H14), you need to modify the mass assigment accordingly...')
                    print('##########')

                elif ID_string == 'C6' and fft == 'CH2_HEX': #Modification made to recognize Ch2 in Hexane
                    FF.atom_masses[fft] = '14.02668'
                    print('##########')
                    print('Mass of CH2 for recognized molecule = ' +str(FF.atom_masses[fft]))
                    print('If you are not trying to include Hexane (C6H14), you need to modify the mass assigment accordingly...')
                    print('##########')


                else:
                    FF.atom_masses[fft] = mass_key[data['element_symbol']]
#                print(FF.atom_masses)

                if 'hybrid' not in FF.pair_data['style'] and style != FF.pair_data['style']:
                    FF.pair_data['style'] = ' '.join(['hybrid', FF.pair_data['style'], style])
                elif 'hybrid' in FF.pair_data['style'] and style in FF.pair_data['style']:
                    pass
                elif 'hybrid' in FF.pair_data['style'] and style not in FF.pair_data['style']:
                    FF.pair_data['style'] += ' ' + style

######### add new bonds
        used_bonds = []
        ty = nbonds

        ### Created for DMF
        if ID_string == 'C3H7N1O1':     #Hard coded for DMF

            new_bond_types[('C_CH3_dmf', 'H_CH3_dmf')]  = ty + 1
            new_bond_types[('C_CH3_dmf', 'N_dmf')]      = ty + 2
            new_bond_types[('C_COH_dmf', 'N_dmf')]      = ty + 3
            new_bond_types[('C_COH_dmf', 'O_dmf')]      = ty + 4
            new_bond_types[('C_COH_dmf', 'H_COH_dmf')]  = ty + 5


            for e0,e1,data in subG.edges(data=True):

                bonds = constants['bonds']['type']   
                #print(constants)
                fft_i = SG.nodes[e0]['force_field_type']
                fft_j = SG.nodes[e1]['force_field_type']
                # make sure the order corresponds to that in the molecule dictionary
                bond = tuple(sorted([fft_i, fft_j]))
        
                try:
                    style = bonds[bond][0]

                    if bond not in used_bonds:

                        FF.bond_data['params'][new_bond_types[bond]] = list(bonds[bond])
                        FF.bond_data['comments'][new_bond_types[bond]] = list(bond)

                        used_bonds.append(bond)

                        #print (new_bond_types[bond], bond)


                    if 'hybrid' not in FF.bond_data['style'] and style != FF.bond_data['style']:
                        FF.bond_data['style'] = ' '.join(['hybrid', FF.bond_data['style'], style])

                    elif 'hybrid' in FF.bond_data['style'] and style in FF.bond_data['style']:
                        pass

                    elif 'hybrid' in FF.bond_data['style'] and style not in FF.bond_data['style']:
                        FF.bond_data['style'] += ' ' + style

                    
                    if new_bond_types[bond] in FF.bond_data['all_bonds']:
                        bond_type = new_bond_types[bond]
                        FF.bond_data['count'] = (FF.bond_data['count'][0] + 1, FF.bond_data['count'][1] + 1)
                        FF.bond_data['all_bonds'][bond_type].append((e0,e1))
                    else:
                        bond_type = new_bond_types[bond]
                        FF.bond_data['count'] = (FF.bond_data['count'][0] + 1, FF.bond_data['count'][1] + 1)
                        FF.bond_data['all_bonds'][bond_type] = [(e0,e1)]

                    # Organizing the bond_data in numeric order 12 = (ty+5) angle types + 1 (python)
                    sorted_bond_params = {key: FF.bond_data['params'][key] for key in range(1, ty+5+1) if key in  FF.bond_data['params']}
                    FF.bond_data['params'].clear()
                    FF.bond_data['params'].update(sorted_bond_params)

                    sorted_bond_comment = {key: FF.bond_data['comments'][key] for key in range(1, ty+5+1) if key in  FF.bond_data['comments']}
                    FF.bond_data['comments'].clear()
                    FF.bond_data['comments'].update(sorted_bond_comment)

                    #print(FF.bond_data['params'])#.keys())
                    #print(FF.bond_data['comments'].keys())
                    #print('')

                except KeyError:
                    pass

        ### This is the general bond type allocation (complex molecules required settting like the above)
        else:

            for e0,e1,data in subG.edges(data=True):

                        #print(constants['bonds'])
                        bonds = constants['bonds']['type']          #This requires the FF to include types. If the types are removed, it works for molecules with 1 bond type.
                        
                        fft_i = SG.nodes[e0]['force_field_type']
                        fft_j = SG.nodes[e1]['force_field_type']
                        # make sure the order corresponds to that in the molecule dictionary
                        bond = tuple(sorted([fft_i, fft_j]))
                        #print(bond)
                
                        try:
                            style = bonds[bond][0]

                            if bond not in used_bonds:

                                ty = ty + 1
                                new_bond_types[bond] = ty
                                FF.bond_data['params'][ty] = list(bonds[bond])
                                FF.bond_data['comments'][ty] = list(bond)

                                used_bonds.append(bond)


                            if 'hybrid' not in FF.bond_data['style'] and style != FF.bond_data['style']:
                                FF.bond_data['style'] = ' '.join(['hybrid', FF.bond_data['style'], style])

                            elif 'hybrid' in FF.bond_data['style'] and style in FF.bond_data['style']:
                                pass

                            elif 'hybrid' in FF.bond_data['style'] and style not in FF.bond_data['style']:
                                FF.bond_data['style'] += ' ' + style

                            #print(FF.bond_data)
                            
                            if ty in FF.bond_data['all_bonds']:
                                bond_type = new_bond_types[bond]
                                FF.bond_data['count'] = (FF.bond_data['count'][0] + 1, FF.bond_data['count'][1] + 1)
                                FF.bond_data['all_bonds'][bond_type].append((e0,e1))
                            else:
                                bond_type = new_bond_types[bond]
                                FF.bond_data['count'] = (FF.bond_data['count'][0] + 1, FF.bond_data['count'][1] + 1)
                                FF.bond_data['all_bonds'][bond_type] = [(e0,e1)]

                            #print(FF.bond_data)

                        except KeyError:
                            pass
        

        
######### add new angles
        used_angles = []
        ty = nangles

        ### Created for DMF
        if ID_string == 'C3H7N1O1':
            #print(ID_string)

            new_angle_types[('H_CH3_dmf', 'C_CH3_dmf', 'H_CH3_dmf')]    = ty + 1
            new_angle_types[('H_CH3_dmf', 'C_CH3_dmf', 'N_dmf')]        = ty + 2
            new_angle_types[('C_CH3_dmf', 'N_dmf', 'C_CH3_dmf')]        = ty + 3
            new_angle_types[('C_CH3_dmf', 'N_dmf', 'C_COH_dmf')]        = ty + 4
            new_angle_types[('N_dmf', 'C_COH_dmf', 'O_dmf')]            = ty + 5
            new_angle_types[('H_COH_dmf', 'C_COH_dmf', 'N_dmf')]        = ty + 6  
            new_angle_types[('H_COH_dmf', 'C_COH_dmf', 'O_dmf')]        = ty + 7

            #print('')
            #print('next molecule')
            #print('')

            for name,data in subG.nodes(data=True):

                angles = constants['angles']
                #print(angles)
                nbors = list(subG.neighbors(name))
                #print(nbors)

                for comb in combinations(nbors, 2):

                    j = name
                    i, k = comb
                    fft_i = subG.nodes[i]['force_field_type']
                    fft_j = subG.nodes[j]['force_field_type']
                    fft_k = subG.nodes[k]['force_field_type']

                    angle = sorted((fft_i, fft_k))
                    angle = (angle[0], fft_j, angle[1])
                    #print(angle) 

                    try:
        
                        style = angles[angle][0]
                        #print(style)
                        FF.angle_data['count'] = (FF.angle_data['count'][0] + 1, FF.angle_data['count'][1])
                        #print(FF.angle_data['count'])
                        
                        if angle not in used_angles:
                            
                            #print(new_angle_types)
                            FF.angle_data['count'] = (FF.angle_data['count'][0], FF.angle_data['count'][1] + 1)
                            FF.angle_data['params'][new_angle_types[angle]] = list(angles[angle])
                            FF.angle_data['comments'][new_angle_types[angle]] = list(angle)

                            used_angles.append(angle)

                        #print(angle, new_angle_types[angle])
        
                        if 'hybrid' not in FF.angle_data['style'] and style != FF.angle_data['style']:
                            FF.angle_data['style'] = ' '.join(['hybrid', FF.angle_data['style'], style])

                        elif 'hybrid' in FF.angle_data['style'] and style in FF.angle_data['style']:
                            pass
                            
                        elif 'hybrid' in FF.angle_data['style'] and style not in FF.angle_data['style']:
                            FF.angle_data['style'] += ' ' + style
                        
                        if new_angle_types[angle] in FF.angle_data['all_angles']:
                            angle_type = new_angle_types[angle]
                            FF.angle_data['all_angles'][angle_type].append((i,j,k))

                        else:
                            angle_type = new_angle_types[angle]
                            FF.angle_data['all_angles'][angle_type] = [(i,j,k)]
                        
                        #print(FF.angle_data['all_angles'])

                        #print(str(ty+7))

                        # Organizing the angle_data in numeric order 16 = (ty+7) angle types + 1 (python)
                        sorted_angle_params = {key: FF.angle_data['params'][key] for key in range(1, ty+7+1) if key in  FF.angle_data['params']}
                        FF.angle_data['params'].clear()
                        FF.angle_data['params'].update(sorted_angle_params)

                        sorted_angle_comment = {key: FF.angle_data['comments'][key] for key in range(1, ty+7+1) if key in  FF.angle_data['comments']}
                        FF.angle_data['comments'].clear()
                        FF.angle_data['comments'].update(sorted_angle_comment)

                        #print(FF.angle_data['params'])#.keys())
                        #print(FF.angle_data['comments'].keys())
                        #print('') 

                    except KeyError:
                        pass

        ### This is the general bond type allocation (complex molecules required settting like the above)
        else:

            for name,data in subG.nodes(data=True):

                angles = constants['angles']   #This solves molecules with one angle like water
    #            print(angles)
                nbors = list(subG.neighbors(name))
    #            print(nbors)

                for comb in combinations(nbors, 2):

                    j = name
                    i, k = comb
                    fft_i = subG.nodes[i]['force_field_type']
                    fft_j = subG.nodes[j]['force_field_type']
                    fft_k = subG.nodes[k]['force_field_type']

                    angle = sorted((fft_i, fft_k))
                    angle = (angle[0], fft_j, angle[1])
    #                print(angle) 

                    try:
        
                        style = angles[angle][0]
    #                    print(style)
                        FF.angle_data['count'] = (FF.angle_data['count'][0] + 1, FF.angle_data['count'][1])


                        if angle not in used_angles:
                            
                            ty = ty + 1
                            new_angle_types[angle] = ty
 #                           print(new_angle_types)
                            FF.angle_data['count'] = (FF.angle_data['count'][0], FF.angle_data['count'][1] + 1)
                            FF.angle_data['params'][ty] = list(angles[angle])
                            FF.angle_data['comments'][ty] = list(angle)

                            used_angles.append(angle)

#                        print(used_angles)
        
                        if 'hybrid' not in FF.angle_data['style'] and style != FF.angle_data['style']:
                            FF.angle_data['style'] = ' '.join(['hybrid', FF.angle_data['style'], style])

                        elif 'hybrid' in FF.angle_data['style'] and style in FF.angle_data['style']:
                            pass
                            
                        elif 'hybrid' in FF.angle_data['style'] and style not in FF.angle_data['style']:
                            FF.angle_data['style'] += ' ' + style
                        
                        #print(ty, 'asking')
                        if ty in FF.angle_data['all_angles']:
                            angle_type = new_angle_types[angle]
                            #print(new_angle_types[angle])
                            FF.angle_data['all_angles'][angle_type].append((i,j,k))
                        else:
                            angle_type = new_angle_types[angle]
                            FF.angle_data['all_angles'][angle_type] = [(i,j,k)]
                        #print(FF.angle_data)

                    except KeyError:
                        pass

########## add new dihedrals

        used_dihedrals = []
        ty = ndihedrals

        ### Created for DMF
        if ID_string == 'C3H7N1O1':

            new_dihedral_types[('C_CH3_dmf', 'C_CH3_dmf', 'N_dmf', 'H_CH3_dmf')]    = ty + 1
            new_dihedral_types[('C_COH_dmf', 'C_CH3_dmf', 'N_dmf', 'H_CH3_dmf')]    = ty + 2
            new_dihedral_types[('C_CH3_dmf', 'C_COH_dmf', 'N_dmf', 'O_dmf')]        = ty + 3
            new_dihedral_types[('C_CH3_dmf', 'C_COH_dmf', 'N_dmf', 'H_COH_dmf')]    = ty + 4


            real_subG = subG
            #   Checikg graph shape
            # This could be enhanced towards molecular structure recognition
            #nx.draw(subG, with_labels=True)
            #plt.show()
            
            dihe_combs = permutations(real_subG.nodes(),4)
            #dihe_combs = combinations(real_subG.nodes(),4)

            #print('')
            #print('next molecule')
            #print('')

            # Filter out permutations where the reverse is also present
            filtered_dihe_comb = []
            seen = set()
            for comb in dihe_combs:
                reverse_comb = tuple(reversed(comb))
                if reverse_comb not in seen:
                    filtered_dihe_comb.append(comb)
                    seen.add(comb)

            for combo in filtered_dihe_comb:
            
                dihedrals = constants['dihedrals']
#                print(dihedrals)

                i , j, k, l = combo
                #print(combo)

                if real_subG.has_edge(i,j) and real_subG.has_edge(j,k) and real_subG.has_edge(k,l) :

                    fft_i = real_subG.nodes[i]['force_field_type']
                    fft_j = real_subG.nodes[j]['force_field_type']
                    fft_k = real_subG.nodes[k]['force_field_type']
                    fft_l = real_subG.nodes[l]['force_field_type']

                    '''For DMF there are two dihedral types in our represetnation
                        X-CH3-N-X and X-C_COH-N-X 
                        This approach ensure that all the possible combinations with middle bond
                        C3_N and C2_N are assigned correctly.'''

                    middle = tuple(sorted((fft_j,fft_k)))
                    extremes = tuple(sorted((fft_i,fft_l)))
                
                    dihedral = tuple(((extremes[0],middle[0],middle[1],extremes[1])))
                    #print(dihedral)
                    #print(i, j, k, l)

                    try:
                        style = dihedrals[dihedral][0]
                                                    #### Count of all dihe   , Count of indexed parameters or dihe
                        FF.dihedral_data['count'] = (FF.dihedral_data['count'][0] + 1, FF.dihedral_data['count'][1])

                        if dihedral not in used_dihedrals:
                        
                            #ty = ty + 1
                            #new_dihedral_types[dihedral] = ty
                            FF.dihedral_data['count'] = (FF.dihedral_data['count'][0], FF.dihedral_data['count'][1] + 1)
                            FF.dihedral_data['params'][new_dihedral_types[dihedral]] = list(dihedrals[dihedral])
                            FF.dihedral_data['comments'][new_dihedral_types[dihedral]] = list(dihedral)

                            used_dihedrals.append(dihedral)

                        #print(dihedral, new_dihedral_types[dihedral])
    
                        if 'hybrid' not in FF.dihedral_data['style'] and style != FF.dihedral_data['style']:
                            FF.dihedral_data['style'] = ' '.join(['hybrid', FF.dihedral_data['style'], style])
                        elif 'hybrid' in FF.dihedral_data['style'] and style in FF.dihedral_data['style']:
                            pass
                        elif 'hybrid' in FF.dihedral_data['style'] and style not in FF.dihedral_data['style']:
                            FF.dihedral_data['style'] += ' ' + style

                        #print(ty, 'asking')

                        if new_dihedral_types[dihedral] in FF.dihedral_data['all_dihedrals']:
                            dihedral_type = new_dihedral_types[dihedral]
                            FF.dihedral_data['all_dihedrals'][dihedral_type].append((i,j,k,l))
                        else:
                            dihedral_type = new_dihedral_types[dihedral]
                            FF.dihedral_data['all_dihedrals'][dihedral_type] = [(i,j,k,l)]

                        #print(FF.dihedral_data['all_dihedrals'])

                        # Organizing the bond_data in numeric order 8 = (ty+4) dihe types + 1 (python)
                        sorted_dihe_params = {key: FF.dihedral_data['params'][key] for key in range(1, ty+4+1) if key in  FF.dihedral_data['params']}
                        FF.dihedral_data['params'].clear()
                        FF.dihedral_data['params'].update(sorted_dihe_params)

                        sorted_dihe_comment = {key: FF.dihedral_data['comments'][key] for key in range(1, ty+4+1) if key in  FF.dihedral_data['comments']}
                        FF.dihedral_data['comments'].clear()
                        FF.dihedral_data['comments'].update(sorted_dihe_comment)

                        #print(FF.bond_data['params'])#.keys())
                        #print(FF.bond_data['comments'].keys())
                        #print('')

                    except KeyError:
                        pass

#######
        ### Created for Hexane
        elif ID_string == 'C6':

            new_dihedral_types[('CH2_HEX', 'CH2_HEX', 'CH2_HEX', 'CH3_HEX')]    = ty + 1
            new_dihedral_types[('CH2_HEX', 'CH2_HEX', 'CH2_HEX', 'CH2_HEX')]    = ty + 2


            real_subG = subG
            #   Checikg graph shape
            # This could be enhanced towards molecular structure recognition
            #nx.draw(subG, with_labels=True)
            #plt.show()
            
            dihe_combs = permutations(real_subG.nodes(),4)

            #print('')
            #print('next molecule')
            #print('')

            # Filter out permutations where the reverse is also present
            filtered_dihe_comb = []
            seen = set()
            for comb in dihe_combs:
                reverse_comb = tuple(reversed(comb))
                if reverse_comb not in seen:
                    filtered_dihe_comb.append(comb)
                    seen.add(comb)

            for combo in filtered_dihe_comb:
            
                dihedrals = constants['dihedrals']
#                print(dihedrals)

                i , j, k, l = combo
#                print(combo)

                if real_subG.has_edge(i,j) and real_subG.has_edge(j,k) and real_subG.has_edge(k,l) :

                    fft_i = real_subG.nodes[i]['force_field_type']
                    fft_j = real_subG.nodes[j]['force_field_type']
                    fft_k = real_subG.nodes[k]['force_field_type']
                    fft_l = real_subG.nodes[l]['force_field_type']

                    #print(fft_i)

                    '''For C6H14 there are two dihedral types in our represetnation
                        CH3-CH2-CH2-CH2 and CH2-CH2-CH2-CH2 
                        This approach ensure that all the possible combinations with middle bond
                        CH2-CH2 are assigned correctly.'''

                    middle = tuple(sorted((fft_j,fft_k)))
                    extremes = tuple(sorted((fft_i,fft_l)))
                
                    dihedral = tuple(((extremes[0],middle[0],middle[1],extremes[1])))
                    #print(dihedral)
                    #print(i, j, k, l)


                    try:

                        style = dihedrals[dihedral][0]
                                                        #### Count of all dihe   , Count of indexed parameters or dihe
                        FF.dihedral_data['count'] = (FF.dihedral_data['count'][0] + 1, FF.dihedral_data['count'][1])
    
                        if dihedral not in used_dihedrals:
                        
                            ty = ty + 1
                            new_dihedral_types[dihedral] = ty
                            FF.dihedral_data['count'] = (FF.dihedral_data['count'][0], FF.dihedral_data['count'][1] + 1)
                            FF.dihedral_data['params'][ty] = list(dihedrals[dihedral])
                            FF.dihedral_data['comments'][ty] = list(dihedral)

                            used_dihedrals.append(dihedral)
                            #print(used_dihedrals)
    
                        if 'hybrid' not in FF.dihedral_data['style'] and style != FF.dihedral_data['style']:
                            FF.dihedral_data['style'] = ' '.join(['hybrid', FF.dihedral_data['style'], style])
                        elif 'hybrid' in FF.dihedral_data['style'] and style in FF.dihedral_data['style']:
                            pass
                        elif 'hybrid' in FF.dihedral_data['style'] and style not in FF.dihedral_data['style']:
                            FF.dihedral_data['style'] += ' ' + style
                    
                        if ty in FF.dihedral_data['all_dihedrals']:
                            
                            FF.dihedral_data['all_dihedrals'][ty].append((i,j,k,l))
                        else:
                            FF.dihedral_data['all_dihedrals'][ty] = [(i,j,k,l)]

                        #print(FF.dihedral_data['params'])

                    except KeyError:
                        pass


#Include improper 
#################

        used_impropers = []
        ty = nimpropers

        ### Created for DMF
        if ID_string == 'C3H7N1O1':
            #print('')
            #print('next molecule')
            #print('')

            for name,data in real_subG.nodes(data=True):

                impropers = constants['impropers']

                nbors = list(subG.neighbors(name))

                for comb in combinations(nbors, 3):

                    i = name
                    j, k, l = comb
        
                    fft_i = real_subG.nodes[i]['force_field_type']
                    fft_j = real_subG.nodes[j]['force_field_type']
                    fft_k = real_subG.nodes[k]['force_field_type']
                    fft_l = real_subG.nodes[l]['force_field_type']

                    improper = tuple(sorted((fft_i, fft_j, fft_k, fft_l)))
                    improper = (improper[0], improper[1], improper[2], improper[3])
                    #print(improper)

                    try:

                        style = impropers[improper][0]
                                                        #### Count of all impro   , Count of indexed parameters of impro
                        FF.improper_data['count'] = (FF.improper_data['count'][0] + 1, FF.improper_data['count'][1])
    
                        if improper not in used_impropers:
                        
                            ty = ty + 1
                            new_improper_types[improper] = ty
                            FF.improper_data['count'] = (FF.improper_data['count'][0], FF.improper_data['count'][1] + 1)
                            FF.improper_data['params'][ty] = list(impropers[improper])
                            FF.improper_data['comments'][ty] = list(improper)

                            used_impropers.append(improper)
                            #print(used_impropers)
    
                        if 'hybrid' not in FF.improper_data['style'] and style != FF.improper_data['style']:
                            FF.improper_data['style'] = ' '.join(['hybrid', FF.improper_data['style'], style])

                        elif 'hybrid' in FF.improper_data['style'] and style in FF.improper_data['style']:
                            pass
                        elif 'hybrid' in FF.improper_data['style'] and style not in FF.improper_data['style']:
                            FF.improper_data['style'] += ' ' + style

                    
                        if ty in FF.improper_data['all_impropers']:
                            FF.improper_data['all_impropers'][ty].append((i,j,k,l))
                        else:
                            FF.improper_data['all_impropers'][ty] = [(i,j,k,l)]

                    except KeyError:
                        pass     


#########
    
    FF.bond_data['count'] = (FF.bond_data['count'][0], len(FF.bond_data['params']))
    FF.angle_data['count'] = (FF.angle_data['count'][0], len(FF.angle_data['params']))
    FF.dihedral_data['count'] = (FF.dihedral_data['count'][0], len(FF.dihedral_data['params']))
    FF.improper_data['count'] = (FF.improper_data['count'][0], len(FF.improper_data['params']))


### This section is needed for inclusion of all the parameters in the inputfile

    #Water
    if 'tip4p' in FF.pair_data['style']:

        for ty,pair in FF.pair_data['comments'].items():
            fft = pair[0]
            if fft == 'O_w':
                FF.pair_data['O_type'] = ty
            if fft == 'H_w':
                FF.pair_data['H_type'] = ty

        for ty, bond in FF.bond_data['comments'].items():
            if sorted(bond) == ['H_w', 'O_w']:
                FF.pair_data['H2O_bond_type'] = ty

        for ty, angle in FF.angle_data['comments'].items():
            if angle == ['H_w', 'O_w', 'H_w']:
                FF.pair_data['H2O_angle_type'] = ty

        if 'long' in FF.pair_data['style']:
            FF.pair_data['M_site_dist'] = 0.1250 #TIP4P Long or (Ewald)
        elif 'cut' in FF.pair_data['style'] and ff_string == 'TIP4P_2005_cutoff':
            FF.pair_data['M_site_dist'] = 0.1546
        elif 'cut' in FF.pair_data['style'] and ff_string == 'TIP4P_cutoff':
            FF.pair_data['M_site_dist'] = 0.1500

##########
#Solving for Methanol and Hexane
    #if 'TraPPE' in FF.pair_data['style']:
    if ff_string  == 'TraPPE':
        if ID_string == 'C1H1O1':

            for ty,pair in FF.pair_data['comments'].items():
                fft = pair[0]
                #print(FF.pair_data)

                if fft == 'C_CH3OH':
                    FF.pair_data['Ch3_type'] = ty
                if fft == 'O_CH3OH':
                    FF.pair_data['O_type'] = ty
                if fft == 'H_CH3OH':
                    FF.pair_data['H_type'] = ty

            for ty, bond in FF.bond_data['comments'].items():
                if sorted(bond) == ['C_CH3OH', 'O_CH3OH']:
                    FF.pair_data['Ch3_OH_bond_type'] = ty

                if sorted(bond) == ['H_CH3OH', 'O_CH3OH']:
                    FF.pair_data['O_H_bond_type'] = ty

            for ty, angle in FF.angle_data['comments'].items():
                if angle == ['C_CH3OH', 'O_CH3OH', 'H_CH3OH']:
                    FF.pair_data['CH3OH_angle_type'] = ty

#### Solving for Hexane
    #print(FF.pair_data)
        #if 'hybrid' in FF.pair_data['style']:
        if ID_string == 'C6':

            for ty,pair in FF.pair_data['comments'].items():
                fft = pair[0]
                #print(fft)

                if fft == 'CH3_HEX':
                    FF.pair_data['Ch3_type'] = ty
                if fft == 'CH2_HEX':
                    FF.pair_data['Ch2_type'] = ty


            for ty, bond in FF.bond_data['comments'].items():

                if sorted(bond) == ['CH2_HEX', 'CH3_HEX']:
                    FF.pair_data['Ch3_Ch2_bond_type'] = ty

                if sorted(bond) == ['CH2_HEX', 'CH2_HEX']:
                    FF.pair_data['Ch2_Ch2_bond_type'] = ty

            for ty, angle in FF.angle_data['comments'].items():

                if angle == ['CH2_HEX', 'CH2_HEX', 'CH3_HEX']:
                    FF.pair_data['Ch2_Ch2_Ch3_angle_type'] = ty
                if angle == ['CH2_HEX', 'CH2_HEX', 'CH2_HEX']:
                    FF.pair_data['Ch2_Ch2_Ch2_angle_type'] = ty


            for ty, dihedral in FF.dihedral_data['comments'].items():
                #print(dihedral)
                if dihedrals == ['Ga_HEX', 'Ga_HEX', 'Pd_HEX']:
                    FF.pair_data['Ch2_Ch2_Ch3_angle_type'] = ty
                if dihedral == ['Ga_HEX', 'Ga_HEX', 'Ga_HEX']:
                    FF.pair_data['Ch2_Ch2_Ch2_angle_type'] = ty

##########
#### Solving for DMF

    if ff_string  == 'UFF':
        
        for ty,pair in FF.pair_data['comments'].items():
            fft = pair[0]

            if fft == 'H_CH3_dmf':
                FF.pair_data['H_CH3_type'] = ty
            if fft == 'C_CH3_dmf':
                FF.pair_data['C_CH3_type'] = ty
            if fft == 'N_dmf':
                FF.pair_data['N_type'] = ty
            if fft == 'C_COH_dmf':
                FF.pair_data['C_COH_type'] = ty
            if fft == 'O_dmf':
                FF.pair_data['O_type'] = ty
            if fft == 'H_COH_dmf':
                FF.pair_data['H_COH_type'] = ty
#        print(FF.pair_data)

        for ty, bond in FF.bond_data['comments'].items():
            if sorted(bond) == ['C_CH3_dmf', 'H_CH3_dmf']:
                FF.pair_data['H_CH3_bond_type'] = ty

            if sorted(bond) == ['C_CH3_dmf', 'N_dmf']:
                FF.pair_data['N_CH3_bond_type'] = ty

            if sorted(bond) == ['C_COH_dmf', 'N_dmf']:
                FF.pair_data['N_COH_bond_type'] = ty

            if sorted(bond) == ['C_COH_dmf', 'O_dmf']:
                FF.pair_data['O_COH_bond_type'] = ty

            if sorted(bond) == ['C_COH_dmf', 'H_COH_dmf']:
                FF.pair_data['H_COH_bond_type'] = ty


        for ty, angle in FF.angle_data['comments'].items():
            if angle == ['H_CH3_dmf', 'C_CH3_dmf', 'N_dmf']:
                FF.pair_data['H_CH3_N_angle_type'] = ty

            if angle == ['H_CH3_dmf', 'C_CH3_dmf', 'H_CH3_dmf']:
                FF.pair_data['H_CH3_H_angle_type'] = ty

            if angle == ['C_CH3_dmf', 'N_dmf', 'C_CH3_dmf']:
                FF.pair_data['CH3_N_CH3__angle_type'] = ty

            if angle == ['C_CH3_dmf', 'N_dmf', 'C_COH_dmf']:
                FF.pair_data['COH_N_CH3_angle_type'] = ty

            if angle == ['N_dmf', 'C_COH_dmf', 'O_dmf']:
                FF.pair_data['N_COH_O_angle_type'] = ty

            if angle == ['H_COH_dmf', 'C_COH_dmf', 'N_dmf']:
                FF.pair_data['H_COH_N_angle_type'] = ty

            if angle == ['H_COH_dmf', 'C_COH_dmf', 'O_dmf']:
                FF.pair_data['H_COH_O_angle_type'] = ty

        #print(FF.dihedral_data['comments'])
        for ty, dihedral in FF.dihedral_data['comments'].items():

            if dihedral == ['C_COH_dmf', 'C_CH3_dmf', 'N_dmf', 'H_CH3_dmf']:
                FF.pair_data['COH_N_CH3_H_dihe_type'] = ty
    
            if dihedral == ['C_CH3_dmf', 'C_CH3_dmf', 'N_dmf', 'H_CH3_dmf']:
                FF.pair_data['CH3_N_CH3_H_dihe_type'] = ty

            if dihedral == ['C_CH3_dmf', 'C_COH_dmf', 'N_dmf', 'O_dmf']:
                FF.pair_data['CH3_N_COH_O_dihe_type'] = ty

            if dihedral == ['C_CH3_dmf', 'C_COH_dmf', 'N_dmf', 'H_COH_dmf']:
                FF.pair_data['CH3_N_COH_H_dihe_type'] = ty

        #print(FF.improper_data['comments'].items())

        for ty, improper in FF.improper_data['comments'].items():
            #print(improper)
            if improper == ['C_COH_dmf', 'H_COH_dmf', 'N_dmf', 'O_dmf']:

                FF.pair_data['COH_H_O_N_impro_type'] = ty

    else:
        #except KeyError:
        pass
##########


def update_potential(potential_data, new_potential_params, potential_coeff):

    write_instyles = False
    add_styles = set([new_potential_params[ty]['style'] for ty in new_potential_params])
    for ABS in add_styles:
        if ABS not in potential_data['style'] and 'hybrid' in potential_data['style']:
            potential_data['style'] = potential_data['style'] + ' ' + ABS
            write_instyles = True
        if ABS not in potential_data['style'] and 'hybrid' not in potential_data['style']:
            potential_data['style'] = 'hybrid ' + potential_data['style'] + ' ' + ABS
            write_instyles = True
        if ABS in potential_data['style'] and 'hybrid' in potential_data['style']:
            potential_data['style'] = potential_data['style'] #+ ' ' + ABS
            write_instyles = True            
        else:
            pass

    if write_instyles:
        instyles = {ty:' ' + new_potential_params[ty]['style'] for ty in new_potential_params}

    else:
        instyles = {ty:'' for ty in new_potential_params}

    potential_data['infile_add_lines'] = []
    for ty,data in new_potential_params.items():
        strparams = ' '.join([str(p) for p in data['params']])
        potential_data['infile_add_lines'].append(potential_coeff + str(ty) + instyles[ty] + ' ' + strparams + ' ' + data['comments'])


def include_molecule_file(FF, maxIDs, add_molecule):

    max_atom_ty, max_bond_ty, max_angle_ty, max_dihedral_ty, max_improper_ty = maxIDs
    molname, model, N = add_molecule

    if molname in ('WATER','water', 'Water', 'H2O', 'h2o'):
        
        molfile, LJ_params, bond_params, angle_params, molnames, mass_dict, M_site_dist, extra_types = WMF.water(max_atom_ty, max_bond_ty, max_angle_ty, model=model)
        dihedral_params = None
        improper_params = None
        FF.pair_data['special_bonds'] = 'lj 0.0 0.0 1.0 coul 0.0 0.0 0.0'  # This special bond is used to be consistent through all the molecules and the MOFs
        FF.pair_data['O_type'] = max_atom_ty + 1
        FF.pair_data['H_type'] = max_atom_ty + 2
        FF.pair_data['H2O_bond_type'] = max_bond_ty + 1
        FF.pair_data['H2O_angle_type'] = max_angle_ty + 1
        FF.pair_data['M_site_dist'] = M_site_dist

## Inclusion of new molecules
#Methanol

    if molname in ('METHANOL','methanol', 'Methanol', 'CH3OH', 'ch3oh'):
        
        molfile, LJ_params, bond_params, angle_params, molnames, mass_dict, M_site_dist, extra_types = WMF.methanol(max_atom_ty, max_bond_ty, max_angle_ty, model=model)
        dihedral_params = None
        improper_params = None
        FF.pair_data['special_bonds'] = 'lj 0.0 0.0 1.0 coul 0.0 0.0 0.0'
        FF.pair_data['Ch3_type'] = max_atom_ty + 1
        FF.pair_data['O_type'] = max_atom_ty + 2
        FF.pair_data['H_type'] = max_atom_ty + 3
        FF.pair_data['Ch3_OH_bond_type'] = max_bond_ty + 1
        FF.pair_data['O_H_bond_type'] = max_bond_ty + 2
        FF.pair_data['Ch3OH_angle_type'] = max_angle_ty + 1

#Hexane
    if molname in ('HEXANE','hexane', 'Hexane', 'C6H14', 'c6h14'):
        
        molfile, LJ_params, bond_params, angle_params, dihedral_params , molnames, mass_dict, M_site_dist, extra_types = WMF.hexane(max_atom_ty, max_bond_ty, max_angle_ty, max_dihedral_ty, model=model)
        improper_params = None
        FF.pair_data['special_bonds'] = 'lj 0.0 0.0 1.0 coul 0.0 0.0 0.0'
        FF.pair_data['Ch3_type'] = max_atom_ty + 1
        FF.pair_data['Ch2_type'] = max_atom_ty + 2
        FF.pair_data['Chx_Chy_bond_type'] = max_bond_ty + 1
        FF.pair_data['Chx_Chy_Chz_angle_type'] = max_angle_ty + 1
        FF.pair_data['Chw_Chx_Chy_Chz_dihe_type'] = max_dihedral_ty + 1

#DMF
    if molname in ('DMF','Dmf','dmf'):
        molfile, LJ_params, bond_params, angle_params, dihedral_params, improper_params, molnames, mass_dict, M_site_dist, extra_types = WMF.dmf(max_atom_ty, max_bond_ty, max_angle_ty, max_dihedral_ty, max_improper_ty, model=model)
        
        #improper_params = None
        FF.pair_data['special_bonds'] = 'lj 0.0 0.0 1.0 coul 0.0 0.0 0.0'
        #   Atom types
        FF.pair_data['C_CH3_type'] = max_atom_ty + 1
        FF.pair_data['C_COH_type'] = max_atom_ty + 2
        FF.pair_data['N_type'] = max_atom_ty + 3
        FF.pair_data['O_type'] = max_atom_ty + 4
        FF.pair_data['H_CH3_type'] = max_atom_ty + 5
        FF.pair_data['H_COH_type'] = max_atom_ty + 6
        # Bond types
        FF.pair_data['H_CH3_bond_type'] = max_bond_ty + 1
        FF.pair_data['N_CH3_bond_type'] = max_bond_ty + 2
        FF.pair_data['N_COH_bond_type'] = max_bond_ty + 3
        FF.pair_data['O_COH_bond_type'] = max_bond_ty + 4
        FF.pair_data['H_COH_bond_type'] = max_bond_ty + 5
        #   Angle types
        FF.pair_data['H_CH3_N_angle_type'] = max_angle_ty + 1
        FF.pair_data['H_CH3_H_angle_type'] = max_angle_ty + 2
        FF.pair_data['CH3_N_CH3__angle_type'] = max_angle_ty + 3
        FF.pair_data['COH_N_CH3_angle_type'] = max_angle_ty + 4
        FF.pair_data['N_COH_O_angle_type'] = max_angle_ty + 5
        FF.pair_data['H_COH_N_angle_type'] = max_angle_ty + 6
        FF.pair_data['H_COH_O_angle_type'] = max_angle_ty + 7
        #   Dihedral types
        FF.pair_data['CH3_N_CH3_H_dihe_type'] = max_dihedral_ty + 1
        FF.pair_data['COH_N_CH3_H_dihe_type'] = max_dihedral_ty + 1
        FF.pair_data['CH3_N_COH_O_dihe_type'] = max_dihedral_ty + 2
        FF.pair_data['CH3_N_COH_H_dihe_type'] = max_dihedral_ty + 2
        #   Improper types
        FF.pair_data['COH_H_O_N_impro_type'] = max_improper_ty + 1


################

    add_LJ_style = LJ_params['style']

    if add_LJ_style not in FF.pair_data['style']:
        FF.pair_data['style'] = FF.pair_data['style'] + ' ' + add_LJ_style
        if 'hybrid' not in FF.pair_data['style']:
            FF.pair_data['style'] = 'hybrid ' + FF.pair_data['style']

    for aty, param in LJ_params.items():
        if aty not in ('style', 'comments'):
            FF.pair_data['params'][aty] = param
            FF.pair_data['comments'][aty] = LJ_params['comments'][aty]

    if bond_params != None:
        update_potential(FF.bond_data, bond_params, 'bond_coeff      ')

    if angle_params != None:
        update_potential(FF.angle_data, angle_params, 'angle_coeff     ')

    if dihedral_params != None:
        update_potential(FF.dihedral_data, dihedral_params, 'dihedral_coeff  ')

    if improper_params != None:
        update_potential(FF.improper_data, improper_params, 'improper_coeff  ')


    infile_add_lines = ['molecule        ' + ' '.join(molnames)]
    for atom in mass_dict:
        infile_add_lines.append('mass            ' + str(atom) + ' ' + str(mass_dict[atom]))
    seed0 = randint(1,10000)
    seed1 = randint(1,10000)

    if N > 0:
        create_line = ' '.join([str(N), str(seed0), 'NULL', 'mol', molnames[0], str(seed1), 'units', 'box'])
        infile_add_lines.append('create_atoms    0 random ' + create_line)

    return molfile, infile_add_lines, extra_types

def read_RASPA_pdb(file):

    with open(file, 'r') as pdb:
        pdb = pdb.read()
        pdb = pdb.split('\n')

    atoms = Atoms()
    for line in pdb:
        s = line.split()
        if len(s) > 0 and s[0] == 'ATOM':
            atoms.append(Atom(s[2], np.array([float(c) for c in s[4:7]])))

    return atoms

def read_small_molecule_file(sm_file, system):

    fm = sm_file.split('.')[-1]
    max_ind = max(system['graph'])
    ind = max_ind + 1

    if fm not in ('pdb', 'xyz', 'cif'):
        raise ValueError('only xyz and RASPA pdb formats are supported for small molecule files')

    if fm == 'pdb':
        print('assuming small molecule file is in RASPA pdb format, if not, too bad...')
        atoms = read_RASPA_pdb(sm_file)
    else:
        atoms = read(sm_file, format=fm)

    atoms.set_cell(system['box'])
    SMG = nx.Graph()

    for atom in atoms:
        
        SMG.add_node(ind, element_symbol=atom.symbol, mol_flag='1', index=ind, force_field_type='', cartesian_position=atom.position, 
                     fractional_position=atom.scaled_position, charge=0.0, replication=np.array([0.0,0.0,0.0]), duplicated_version_of=None)

        ind += 1

    system['SM_graph'] = nx.compose(SMG, system['SM_graph']) # don't want to overwrite extra framework species already in the cif
