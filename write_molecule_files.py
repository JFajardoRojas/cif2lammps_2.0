from textwrap import dedent

def water(last_atom_ID, last_bond_ID, last_angle_ID, model='TIP4P_cutoff'):

    ID_O = last_atom_ID + 1
    ID_H = ID_O + 1

    BT = last_bond_ID + 1
    AT = last_angle_ID + 1

    charge_dict = {
    'TIP4P_cutoff': (-1.04000, 0.52000),
    'TIP4P_2005':   (-1.11280, 0.55640),
    'TIP4P_long':   (-1.04844, 0.52422), # This is the default for adding molecules for FE solvation
    'TIP3P_long':   (-0.83000, 0.41500)
    }

    M_site_dist_dict = {
    'TIP4P_cutoff': 0.1500,
    'TIP4P_2005':   0.1546,
    'TIP4P_long':   0.1250, # This is the default for adding molecules for FE solvation
    'TIP3P_long':   None
    }

    LJ_dict = {
    # LAMMPS has a special TIP4P pair_style that automatically adds the M site
    'TIP4P_cutoff': {ID_O: ('lj/cut/tip4p/cut' , 0.15500, 3.15360), ID_H: ('lj/cut/tip4p/cut' , 0.0, 1.0), 'style': 'lj/cut/tip4p/cut' , 'comments': {ID_O:['O_water', 'O_water'], ID_H:['H_water', 'H_water']}},
    'TIP4P_2005':   {ID_O: ('lj/cut/tip4p/long', 0.18520, 3.15890), ID_H: ('lj/cut/tip4p/long', 0.0, 1.0), 'style': 'lj/cut/tip4p/long', 'comments': {ID_O:['O_water', 'O_water'], ID_H:['H_water', 'H_water']}},
    'TIP4P_long':   {ID_O: ('lj/cut/tip4p/long', 0.16275, 3.16435), ID_H: ('lj/cut/tip4p/long', 0.0, 1.0), 'style': 'lj/cut/tip4p/long', 'comments': {ID_O:['O_water', 'O_water'], ID_H:['H_water', 'H_water']}},
    'TIP3P_long':   {ID_O: ('lj/cut/coul/long' , 0.10200, 3.18800), ID_H: ('lj/cut/coul/long' , 0.0, 1.0), 'style': 'lj/cut/coul/long' , 'comments': {ID_O:['O_water', 'O_water'], ID_H:['H_water', 'H_water']}}
    }

    bond_dict = {
    # TIP4P is a rigid model (use fix shake), force constants should just be ~reasonable values
    # TIP3P has force constants if a flexible model is desired
    'TIP4P_cutoff': {BT: {'style':'harmonic', 'params':(400.0, 0.9572), 'comments':'# O_water H_water'}},
    'TIP4P_2005':   {BT: {'style':'harmonic', 'params':(400.0, 0.9572), 'comments':'# O_water H_water'}},
    'TIP4P_long':   {BT: {'style':'harmonic', 'params':(400.0, 0.9572), 'comments':'# O_water H_water'}},
    'TIP3P_long':   {BT: {'style':'harmonic', 'params':(400.0, 0.9572), 'comments':'# O_water H_water'}} 
    }

    angle_dict = {
    # TIP4P is a rigid model (use fix shake), force constants should just be reasonable values
    # TIP3P has force constants if a flexible model is desired
    'TIP4P_cutoff': {AT: {'style':'harmonic', 'params':(50.0, 104.52), 'comments':'# H_water O_water H_water'}},
    'TIP4P_2005':   {AT: {'style':'harmonic', 'params':(50.0, 104.52), 'comments':'# H_water O_water H_water'}},
    'TIP4P_long':   {AT: {'style':'harmonic', 'params':(50.0, 104.52), 'comments':'# H_water O_water H_water'}},
    'TIP3P_long':   {AT: {'style':'harmonic', 'params':(55.0, 104.52), 'comments':'# H_water O_water H_water'}}
    }

    qO,qH = charge_dict[model]
    LJ_params = LJ_dict[model]
    bond_params = bond_dict[model]
    angle_params = angle_dict[model]

    if 'TIP4P' in model:

        molfile = dedent("""# Water molecule. useable for TIP3P or TIP4P in LAMMPS.

3 atoms
2 bonds
1 angles

Coords

1    1.12456   0.09298   1.27452
2    1.53683   0.75606   1.89928
3    0.49482   0.56390   0.65678

Types

1    {ID_O}
2    {ID_H}
3    {ID_H}

Charges

1    {qO}
2    {qH}
3    {qH}

Bonds

1    {BT} 1 2
2    {BT} 1 3

Angles

1    {AT} 2 1 3

Shake Flags

1 1
2 1
3 1

Shake Atoms

1 1 2 3
2 1 2 3
3 1 2 3

Shake Bond Types

1 {BT} {BT} {AT}
2 {BT} {BT} {AT}
3 {BT} {BT} {AT}""".format(**locals())).strip()

    if 'TIP3P' in model:
        
        molfile = dedent("""# Water molecule. useable for TIP3P or TIP4P in LAMMPS.

3 atoms
2 bonds
1 angles

Coords

1    1.12456   0.09298   1.27452
2    1.53683   0.75606   1.89928
3    0.49482   0.56390   0.65678

Types

1    {ID_O}
2    {ID_H}
3    {ID_H}

Charges

1    {qO}
2    {qH}
3    {qH}

Bonds

1    {BT} 1 2
2    {BT} 1 3

Angles

1    {AT} 2 1 3""".format(**locals())).strip()

    mass_dict = {ID_O:15.9994, ID_H:1.00784}
    molnames = ('H2O_mol', 'H2O.txt')

    extra_types = (2,1,1,None,None)

    return molfile, LJ_params, bond_params, angle_params, molnames, mass_dict, M_site_dist_dict[model], extra_types

######

def methanol(last_atom_ID, last_bond_ID, last_angle_ID, model='TraPPE'): #In models different than TraPPE then used impropers and other details required.

    ###
    ID_Ch3 = last_atom_ID + 1
    ID_O = ID_Ch3 + 1
    ID_H = ID_O + 1

    BT_1 = last_bond_ID + 1
    BT_2 = BT_1 + 1

    AT = last_angle_ID + 1

    charge_dict = {
    'TraPPE': (0.265, -0.7, 0.435)
    }

    M_site_dist_dict = {
    'TraPPE': None
    }

    LJ_dict = {
    
    'TraPPE': {ID_Ch3: ('lj/cut/coul/long' , 0.194746, 3.75), ID_O: ('lj/cut/coul/long' , 0.184810, 3.02), ID_H: ('lj/cut/coul/long' , 0.0, 0.0), 'style': 'lj/cut/coul/long' , 'comments': {ID_Ch3:['Ch3_CH3OH','Ch3_CH3OH'],ID_O:['O_CH3OH', 'O_CH3OH'], ID_H:['H_CH3OH', 'H_CH3OH']}}
    }

    bond_dict = {
    # TraPPE is a rigid model (use fix shake), force constants should just be reasonable values
    'TraPPE': {BT_1: {'style':'harmonic', 'params':(400.0, 1.43), 'comments':'# Ch3_CH3OH O_CH3OH'}, BT_2: {'style':'harmonic', 'params':(400.0, 0.945), 'comments':'# O_CH3OH H_CH3OH'}} 
    }

    angle_dict = {
    # TraPPE is a rigid model (use fix shake), force constants should just be reasonable values
    'TraPPE': {AT: {'style':'harmonic', 'params':(55.0, 108.50), 'comments':'# Ch3_CH3OH O_CH3OH H_CH3OH'}}
    }

    qCh3,qO,qH = charge_dict[model]
    LJ_params = LJ_dict[model]
    bond_params = bond_dict[model]
    angle_params = angle_dict[model]

    if 'TraPPE' in model:

        molfile = dedent("""# Methanol molecule. Useable for TraPPE in LAMMPS.

3 atoms
2 bonds
1 angles

Coords

1    1.00000   1.00000   0.00000
2    -0.4160   1.00000   0.00000
3    -0.7010   0.19196   -0.45956

Types

1    {ID_Ch3}
2    {ID_O}
3    {ID_H}

Charges

1    {qCh3}
2    {qO}
3    {qH}

Bonds

1    {BT_1} 1 2
2    {BT_2} 2 3

Angles

1    {AT} 1 2 3

Shake Flags

1 3
2 3
3 3

Shake Atoms

1 2 1 3
2 2 1 3
3 2 1 3

Shake Bond Types

1 {BT_1} {BT_2}
2 {BT_1} {BT_2}
3 {BT_1} {BT_2}""".format(**locals())).strip()

    mass_dict = {ID_Ch3:15.03422, ID_O:15.9994, ID_H:1.00784}
    molnames = ('CH3OH_mol', 'CH3OH.txt')

    extra_types = (3,2,1,None,None)

    return molfile, LJ_params, bond_params, angle_params, molnames, mass_dict, M_site_dist_dict[model], extra_types

############

def hexane(last_atom_ID, last_bond_ID, last_angle_ID, last_dihe_ID, model='TraPPE'): #In models different than TraPPE then used impropers and other details required.

    ###
    ID_Ch3 = last_atom_ID + 1
    ID_Ch2 = ID_Ch3 + 1

    BT = last_bond_ID + 1

    AT = last_angle_ID + 1

    DT = last_dihe_ID + 1

    charge_dict = {
    'TraPPE': (0.0, 0.0)
    }

    M_site_dist_dict = {
    'TraPPE': None
    }

    LJ_dict = {
    # I know that hexane does not have charges (specified as 0.0). I am using this pair style to facilitate the recognition of hybrid pair style useful for the FE_sol
    # calculation. Due to charge 0.0, this change has negligible effect.
    'TraPPE': {ID_Ch3: ('lj/cut' , 0.19461, 3.75), ID_Ch2: ('lj/cut' , 0.09135025, 3.95), 'style': 'lj/cut' , 'comments': {ID_Ch3:['Ch3_C6H14','Ch3_C6H14'],ID_Ch2:['Ch2_C6H14','Ch2_C6H14']}}
    }

    bond_dict = {
    
    # TraPPE is a rigid model (use fix shake), force constants should just be reasonable values
    'TraPPE': {BT: {'style':'harmonic', 'params':(400.0, 1.54), 'comments':'# Ch3_C6H14 Ch2_C6H14'}} 
    }

    angle_dict = {
    'TraPPE': {AT: {'style':'harmonic', 'params':(62.5, 114.0), 'comments':'# Chx_C6H14 Chy_C6H14 Chz_C6H14'}}
    }

    dihe_dict = {
    'TraPPE': {DT: {'style':'opls', 'params':(1.41, -0.27, 3.14, 0.0), 'comments':'# Chw_C6H14 Chx_C6H14 Chy_C6H14 Chz_C6H14'}}
    }

    qCh3,qCh2= charge_dict[model]
    LJ_params = LJ_dict[model]
    bond_params = bond_dict[model]
    angle_params = angle_dict[model]
    dihedral_params = dihe_dict[model]

    if 'TraPPE' in model:

        molfile = dedent("""# Hexane molecule. Useable for TraPPE in LAMMPS. Using Ch3 and Ch2

6 atoms
5 bonds
4 angles
3 dihedrals

Coords

1    1.000  1.00000  0.00000
2    -0.520  1.00000  0.00000
3    -1.079  1.00000  1.42293
4    -2.608  0.99820  1.41965
5    -3.167  1.00000  2.84250
6    -4.687  0.99642  2.84251

Types

1    {ID_Ch3}
2    {ID_Ch2}
3    {ID_Ch2}
4    {ID_Ch2}
5    {ID_Ch2}
6    {ID_Ch3}


Charges

1    {qCh3}
2    {qCh2}
3    {qCh2}
4    {qCh2}
5    {qCh2}
6    {qCh3}

Bonds

1    {BT} 1 2
2    {BT} 2 3
3    {BT} 3 4
4    {BT} 4 5
5    {BT} 5 6

Angles

1    {AT} 1 2 3
2    {AT} 2 3 4
3    {AT} 3 4 5
4    {AT} 4 5 6

Dihedrals

1    {DT} 1 2 3 4
2    {DT} 2 3 4 5
3    {DT} 3 4 5 6

Shake Flags

1 2
2 2
3 2
4 2
5 2
6 2

Shake Atoms

1 1 2 
2 2 3
3 3 4
4 4 5
5 5 6
6 5 6

Shake Bond Types

1 {BT}
2 {BT}
3 {BT}
3 {BT}
3 {BT}
3 {BT}""".format(**locals())).strip()

    mass_dict = {ID_Ch3:15.03422, ID_Ch2:14.0269}
    molnames = ('C6H14_mol', 'C6H14.txt')

    extra_types = (2,1,1,1,None)

    return molfile, LJ_params, bond_params, angle_params, dihedral_params, molnames, mass_dict, M_site_dist_dict[model], extra_types



############

def dmf(last_atom_ID, last_bond_ID, last_angle_ID, last_dihe_ID, last_improp_ID, model='UFF'): 

    ### This molecule is modeled using UFF
    ID_Ch3 = last_atom_ID + 1
    ID_Ch2 = ID_Ch3 + 1
    ID_N = ID_Ch2 + 1
    ID_O = ID_N + 1
    ID_H = ID_O + 1
    ID_H1 = ID_H + 1

    BT_1 = last_bond_ID + 1
    BT_2 = BT_1 + 1
    BT_3 = BT_2 + 1
    BT_4 = BT_3 + 1
    BT_5 = BT_4 + 1

    AT_1 = last_angle_ID + 1
    AT_2 = AT_1 + 1
    AT_3 = AT_2 + 1
    AT_4 = AT_3 + 1
    AT_5 = AT_4 + 1
    AT_6 = AT_5 + 1
    AT_7 = AT_6 + 1

    DT_1 = last_dihe_ID + 1
    DT_2 = DT_1 + 1

    IMP = last_improp_ID + 1


    charge_dict = {
    # The order in chargers is Ch3, Ch2, N, O, H, H1 
    'UFF': (0.1429, 0.5924, -1.1192, -0.4020, 0.0866, 0.1234) #0.1260 Original but a difference of 0.0026 is observed
    }


    M_site_dist_dict = {
    'UFF': None
    }

    LJ_dict = {

    'UFF': {
    ID_Ch3: ('lj/cut/coul/long', 0.10500, 3.43085), # This is in reality C connected to 3 hydorgens. But, it is not a pseudo-atom.
    ID_Ch2: ('lj/cut/coul/long', 0.10500, 3.43085), # This is in reality C connected to H and O rather than Ch2, but I preserve the nomenclature to avoid multiple changes in cif2lammps
    ID_N: ('lj/cut/coul/long', 0.06900, 3.26069),
    ID_O: ('lj/cut/coul/long', 0.06000, 3.11815),
    ID_H: ('lj/cut/coul/long', 0.04400, 2.57113),
    ID_H1: ('lj/cut/coul/long', 0.04400, 2.57113),
    'style': 'lj/cut/coul/long',
    'comments': {
    ID_Ch3:['Ch3_dmf','Ch3_dmf'],
    ID_Ch2:['Ch2_dmf','Ch2_dmf'],
    ID_N:['N_dmf','N_dmf'],
    ID_O:['O_dmf','O_dmf'],
    ID_H:['H_dmf','H_dmf'],
    ID_H1:['H1_dmf','H1_dmf']}
            }
    }

    bond_dict = {
    'UFF': {
    BT_1: {'style':'harmonic', 'params':(331.06939, 1.10940), 'comments':'# Ch3_dmf H_dmf'},
    BT_2: {'style':'harmonic', 'params':(528.59298, 1.45107), 'comments':'# Ch3_dmf N_dmf'},
    BT_3: {'style':'harmonic', 'params':(556.75624, 1.42618), 'comments':'# Ch2_dmf N_dmf'},
    BT_4: {'style':'harmonic', 'params':(796.39449, 1.22393), 'comments':'# Ch2_dmf O_dmf'},
    BT_5: {'style':'harmonic', 'params':(354.48390, 1.08442), 'comments':'# Ch2_dmf H1_dmf'}} 
    }

    angle_dict = {
    'UFF': {
    AT_1: {'style':'fourier', 'params':(169.76232, 0.34375, 0.37500, 0.28125), 'comments':'# H_dmf Ch3_dmf N_dmf'},
    AT_2: {'style':'fourier', 'params':( 75.49653, 0.34375, 0.37500, 0.28125), 'comments':'# H_dmf Ch3_dmf H_dmf'},
    AT_3: {'style':'fourier', 'params':(260.86939, 0.31751, 0.31322, 0.27250), 'comments':'# Ch3_dmf N_dmf Ch3_dmf'},
    AT_4: {'style':'fourier', 'params':(267.65751, 0.31751, 0.31322, 0.27250), 'comments':'# Ch2_dmf N_dmf Ch3_dmf'},
    AT_5: {'style':'cosine/periodic', 'params':(199.35631, -1, 3), 'comments':'# N_dmf Ch2_dmf O_dmf'},
    AT_6: {'style':'cosine/periodic', 'params':( 71.40186, -1, 3), 'comments':'# H1_dmf Ch2_dmf N_dmf'},
    AT_7: {'style':'cosine/periodic', 'params':( 84.67237, -1, 3), 'comments':'# H1_dmf Ch2_dmf O_dmf'}}
    }

    dihe_dict = {
    'UFF': {
    DT_1: {'style':'harmonic', 'params':(0.08138,  1, 3), 'comments':'# X Ch3_dmf N_dmf X'},
    DT_2: {'style':'harmonic', 'params':(0.25000, -1, 6), 'comments':'# X Ch2_dmf N_dmf X'}}
    }

    impro_dict = {
    'UFF': {
    IMP: {'style':'fourier', 'params':(2.0000, 1.0000, -1.0000, 0.0000), 'comments':'# XX Ch3_dmf N_dmf XX'}}
    }

    qCh3, qCh2, qN, qO, qH, qH1 = charge_dict[model]
    LJ_params = LJ_dict[model]
    bond_params = bond_dict[model]
    angle_params = angle_dict[model]
    dihedral_params = dihe_dict[model]
    improper_params = impro_dict[model]

    if 'UFF' in model:

        molfile = dedent("""# DMF molecule. Useable for UFF in LAMMPS.

12 atoms
11 bonds
18 angles
16 dihedrals
1 impropers

Coords

1    4.0000  4.0000  3.0000         
2    2.8840  4.0000  3.0000
3    2.1620  4.0000  4.3040
4    2.5370  3.1370  2.3900
5    2.5260  4.8620  2.3910
6    0.6770  4.0010  4.3040 
7    2.8310  4.1760  5.4350
8    0.1380  4.0000  5.2800
9    0.3200  4.8650  3.6960
10   0.3190  3.1400  3.6940
11   3.3250  3.2100  6.0060
12   2.9810  5.1590  5.9140

Types

1    {ID_H}
2    {ID_Ch3}
3    {ID_N}
4    {ID_H}
5    {ID_H}
6    {ID_Ch3}
7    {ID_Ch2}
8    {ID_H}
9    {ID_H}
10   {ID_H}
11   {ID_O}
12   {ID_H1}

Charges

1    {qH}
2    {qCh3}
3    {qN}
4    {qH}
5    {qH}
6    {qCh3}
7    {qCh2}
8    {qH}
9    {qH}
10   {qH}
11   {qO}
12   {qH1}

Bonds

1    {BT_1} 1 2
2    {BT_1} 2 4
3    {BT_1} 2 5
4    {BT_2} 2 3
5    {BT_2} 6 3
6    {BT_1} 6 8
7    {BT_1} 6 10
8    {BT_1} 6 9
9    {BT_3} 7 3
10   {BT_4} 7 11
11   {BT_5} 7 12

Angles

1    {AT_1} 1 2 3
2    {AT_2} 1 2 4
3    {AT_2} 1 2 5
4    {AT_3} 2 3 6
5    {AT_4} 2 3 7
6    {AT_1} 3 6 8
7    {AT_1} 3 6 9
8    {AT_1} 3 6 10
9    {AT_5} 3 7 11
10   {AT_6} 3 7 12
11   {AT_4} 6 3 7
12   {AT_1} 3 2 4
13   {AT_2} 9 6 10
14   {AT_1} 3 2 5
15   {AT_2} 4 2 5
16   {AT_2} 8 6 10
17   {AT_2} 8 6 9
18   {AT_7} 11 7 12

Dihedrals

1    {DT_1} 6 3 2 1
2    {DT_1} 8 6 3 2
3    {DT_2} 11 7 3 2
4    {DT_1} 6 3 2 4
5    {DT_1} 10 6 3 7
6    {DT_1} 7 3 2 5
7    {DT_1} 7 3 2 1
8    {DT_2} 12 7 3 6
9    {DT_1} 6 3 2 5
10   {DT_1} 7 3 2 4
11   {DT_1} 9 6 3 2
12   {DT_2} 12 7 3 2
13   {DT_1} 8 6 3 7
14   {DT_1} 10 6 3 2
15   {DT_1} 9 6 3 7
16   {DT_2} 11 7 3 6

Impropers

1    {IMP} 7 3 11 12""".format(**locals())).strip()

    mass_dict = {ID_Ch3:12.0107, ID_Ch2:12.0107,ID_N:14.0067,ID_O:15.9994,ID_H:1.00784,ID_H1:1.00784}
    molnames = ('DMF_mol', 'DMF.txt')

    extra_types = (6,5,7,2,1)

    return molfile, LJ_params, bond_params, angle_params, dihedral_params, improper_params, molnames, mass_dict, M_site_dist_dict[model], extra_types



#### More molecules could be integrated.