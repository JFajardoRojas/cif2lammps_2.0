TraPPE =  {
    'O2': {
        'pair': {'style': 'lj/cut/coul/long', 'vdW': {'O_O2': (0.0974,3.02), 'O_com': (0.0,0.0)}, 'charges': {'O_O2': -0.113, 'O_com': 0.226}},
        'bonds': {'type':{('O_O2', 'O_com'): ('harmonic',100.0,0.604)}}, # molecule should be kept rigid, force constants don't matter
        'angles': {('O_O2', 'O_com', 'O_O2'): ('harmonic',100.0,180.0)}, # molecule should be kept rigid, force constants don't matter
        'dihedrals': None,
        'impropers': None
    },
    'N2': {
        'pair': {},
        'bonds': {},
        'angles': {},
        'dihedrals': None,
        'impropers': None
    },
    'H2O1': {
        'pair': {},
        'bonds': {},
        'angles': {},
        'dihedrals': None,
        'impropers': None
    },    
    'C1H1O1': {
        'pair': {'style': 'lj/cut/coul/cut', 'vdW': {'C_CH3OH': (0.19461575, 3.75), 'H_CH3OH': (0.0, 1.0), 'O_CH3OH': (0.184686375, 3.02)}, 'charges': {'C_CH3OH': 0.265, 'H_CH3OH': 0.435, 'O_CH3OH': -0.7}}, 
        'bonds': {'type':{('C_CH3OH', 'O_CH3OH'): ('harmonic',400.0, 1.43), ('H_CH3OH', 'O_CH3OH'): ('harmonic',400.0, 0.945)}}, # molecule should be kept rigid, force constants don't matter
        'angles': {('C_CH3OH', 'O_CH3OH', 'H_CH3OH'): ('harmonic',55.0, 108.50)}, # molecule should be kept rigid, force constants don't matter
        'dihedrals': None,
        'impropers': None
    },    
    'C6': { #To avoid the complications of ASE Ch3 is represented by Pd and Ch2 is represented by Ga.
        'pair': {'style': 'lj/cut/coul/cut', 'vdW': {'CH3_HEX': (0.19461, 3.75), 'CH2_HEX': (0.09135025, 3.95)}, 'charges': {'CH3_HEX': 0.0, 'CH2_HEX': 0.0}},
        'bonds': {'type':{('CH2_HEX', 'CH3_HEX'): ('harmonic',400.0, 1.54),('CH2_HEX', 'CH2_HEX'): ('harmonic',400.0, 1.54)}}, # molecule should be kept rigid, force constants don't matter
        'angles': {('CH2_HEX', 'CH2_HEX', 'CH3_HEX'): ('harmonic',62.5, 114),('CH2_HEX', 'CH2_HEX', 'CH2_HEX'): ('harmonic',62.5, 114)},
        'dihedrals': {('CH2_HEX', 'CH2_HEX', 'CH2_HEX','CH3_HEX'): ('opls',1.41, -0.27,3.14,0.0),('CH2_HEX', 'CH2_HEX', 'CH2_HEX','CH2_HEX'): ('opls',1.41, -0.27,3.14,0.0)},
        'impropers': None
    }
}
# It is needed to include the soft potential, but then it is necessary to include a module to select the applications
# Soft potential does not make sense if the molecule is not created or faded.
UFF =  {
    'C3H7N1O1': {        # This is the string representation of DMF
        'pair': {'style': 'lj/cut/coul/long', 'vdW': {'C_CH3_dmf': (0.10500, 3.43085), 'C_COH_dmf': (0.10500, 3.43085), 'N_dmf': (0.06900, 3.26069), 'O_dmf':(0.06, 3.11815), 'H_CH3_dmf':(0.044, 2.57113), 'H_COH_dmf':(0.044, 2.57113)}, 
        'charges': {'C_CH3_dmf': 0.1429, 'C_COH_dmf': 0.5924, 'N_dmf': -1.1192, 'O_dmf': -0.4020, 'H_CH3_dmf': 0.0866, 'H_COH_dmf': 0.1234}},
        'bonds': {'type':
        {
        ('C_CH3_dmf','H_CH3_dmf'):  ('harmonic', 331.06939, 1.10940), 
        ('C_CH3_dmf','N_dmf'):      ('harmonic', 528.59298, 1.45107),
        ('C_COH_dmf','N_dmf'):      ('harmonic', 556.75624, 1.42618),
        ('C_COH_dmf','H_COH_dmf'):  ('harmonic', 354.48390, 1.08442),
        ('C_COH_dmf','O_dmf'):      ('harmonic', 796.39449, 1.22393) }},
        
        'angles': {
        ('H_CH3_dmf', 'C_CH3_dmf', 'H_CH3_dmf'):    ('fourier', 75.49653, 0.34375, 0.37500, 0.28125),
        ('H_CH3_dmf', 'C_CH3_dmf', 'N_dmf'):        ('fourier', 169.76232, 0.34375, 0.37500, 0.28125),
        ('C_CH3_dmf', 'N_dmf', 'C_CH3_dmf'):        ('fourier', 260.86939, 0.31751, 0.31322, 0.27250),
        ('C_CH3_dmf', 'N_dmf', 'C_COH_dmf'):        ('fourier', 267.65751, 0.31751, 0.31322, 0.27250),
        ('H_COH_dmf', 'C_COH_dmf', 'N_dmf'):        ('cosine/periodic', 71.40186, -1, 3),
        ('N_dmf', 'C_COH_dmf', 'O_dmf'):            ('cosine/periodic', 199.35631, -1, 3),
        ('H_COH_dmf', 'C_COH_dmf', 'O_dmf'):        ('cosine/periodic', 84.67237, -1, 3)}, 
        
        'dihedrals': {  #Due to the recognition protocol. The 16 Dihedrals need to be allocated but there are only
                        # Two set of parameters:
        ('C_COH_dmf', 'C_CH3_dmf', 'N_dmf', 'H_CH3_dmf'):   ('harmonic', 0.08138, 1, 3),  #The middle pair is writen alphabetically for coding purposes
        ('C_CH3_dmf', 'C_CH3_dmf', 'N_dmf', 'H_CH3_dmf'):   ('harmonic', 0.08138, 1, 3),
        ('C_CH3_dmf', 'C_COH_dmf', 'N_dmf', 'H_COH_dmf'):   ('harmonic', 0.25000, -1, 6),  #The middle pair is writen alphabetically for coding purposes
        ('C_CH3_dmf', 'C_COH_dmf', 'N_dmf', 'O_dmf'):       ('harmonic', 0.25000, -1, 6)},
        
        'impropers': {
        ('C_COH_dmf', 'H_COH_dmf', 'N_dmf', 'O_dmf'): ('fourier', 2.00000, 1.00000, -1.00000, 0.00000, 1)},
    }
}


TIP4P_2005_long =  {
# this is TIP4P/2005 water, should be used with long-range electrostatics with 8.5 Å cutoff and fix/shake
# keep in mind that using any long pair_style in lammps will include long-range electrostatics FOR ALL ATOMS in the simulation
    'H2O1': {
        'pair': {'style': 'lj/cut/tip4p/long', 'vdW': {'H_w': (0.0,1.0), 'O_w': (0.162750, 3.164350)}, 'charges': {'H_w': 0.52422, 'O_w': -1.04844}},
        'bonds': {'type':{('H_w', 'O_w'): ('harmonic', 400.0, 0.9572)}},
        'angles': {('H_w', 'O_w', 'H_w'): ('harmonic', 50.0, 104.52)},
        'dihedrals': None,
        'impropers': None
    },
    'Cl1': {
        'pair': {'style': 'lj/cut/tip4p/long', 'vdW': {'Cl_Cl1': (0.22700, 3.51638)}, 'charges': {'Cl_Cl1': -1.0}},
        'bonds': None,
        'angles': None,
        'dihedrals': None,
        'impropers': None
    }
}

TIP4P_2005_cutoff =  {
# this is TIP4P/2005 water but with no long range electrostatics | This is not completely implemented, requires to includ data in write_molecule
    'H2O1': {
        'pair': {'style': 'lj/cut/tip4p/cut', 'vdW': {'H_w': (0.0,0.0), 'O_w': (0.1852, 3.1589)}, 'charges': {'H_w': 0.5564, 'O_w': -1.1128}},
        'bonds': {('H_w', 'O_w'): ('harmonic', 450.0, 0.9572)},
        'angles': {('H_w', 'O_w', 'H_w'): ('harmonic', 55.0, 104.52)},
        'dihedrals': None,
        'impropers': None
    },
    'Cl1': {
        'pair': {'style': 'lj/cut/tip4p/long', 'vdW': {'Cl_Cl1': (0.22700, 3.51638)}, 'charges': {'Cl_Cl1': -1.0}},
        'bonds': None,
        'angles': None,
        'dihedrals': None,
        'impropers': None
    }
}

TIP4P_cutoff =  {
# this is the original TIP4P water model
    'H2O1': {
        'pair': {'style': 'lj/cut/tip4p/cut', 'vdW': {'H_w': (0.0,0.0), 'O_w': (0.1550, 3.1536)}, 'charges': {'H_w': 0.5200, 'O_w': -1.040}},
        'bonds': {('H_w', 'O_w'): ('harmonic', 450.0, 0.9572)},
        'angles': {('H_w', 'O_w', 'H_w'): ('harmonic', 55.0, 104.52)},
        'dihedrals': None,
        'impropers': None
    },
    'Cl1': {
        'pair': {'style': 'lj/cut/tip4p/cut', 'vdW': {'Cl_Cl1': (0.22700, 3.51638)}, 'charges': {'Cl_Cl1': -1.0}},
        'bonds': None,
        'angles': None,
        'dihedrals': None,
        'impropers': None
    }
}

Ions =  {
    'Cl1': {
        'pair': {'style': 'lj/cut/coul/long', 'vdW': {'Cl_Cl1': (0.22700, 3.51638)}, 'charges': {'Cl_Cl1': -1.0}},
        'bonds': None,
        'angles': None,
        'dihedrals': None,
        'impropers': None
    }
}