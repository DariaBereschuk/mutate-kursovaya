import sys
import math
import numpy as np

#РАССЧЕТ УГЛОВ БЛОК 1


def read_pdb(filepath):
    with open(filepath,"r") as f:
        lines=f.readlines()
        atoms = [line.rstrip('\n') for line in lines if line.startswith(('ATOM', 'HETATM'))]
        return atoms

def atoms_in_residue(atoms, chain, res_num):
     res_atoms=[]
     for line in atoms:
         if len(line)<26:
             continue
         line_chain=line[21:22]
         line_res_num=line[22:26].strip()
         if line_chain==chain and res_num==line_res_num:
             res_atoms.append(line)
     return res_atoms

def atoms_xyz(res_atoms, atom_name):
    for line in res_atoms:
        name=line[12:16].strip()
        if name==atom_name:
            x=float(line[30:38])
            y=float(line[38:46])
            z=float(line[46:54])
            return np.array([x,y,z])
        
    return None

def dihedral(p0, p1, p2, p3):
    b0=p0-p1 #вектор связи p1-p0
    b1=p2-p1 #вектор оси вращения, общая связь двух плоскостей
    b2=p3-p2 #вектор связи p2-p3, лежит во второй плоскости
    b1/=np.linalg.norm(b1) #только направление
    v=b0-np.dot(b0, b1)*b1 #проекиц я вектора б0 на плоскость перпендикулярную оси
    w=b2-np.dot(b2,b1)*b1 #вектора б2
    x=np.dot(v, w) #кос угла между проекциями 
    y=np.dot(np.cross(b1,v),w)#синус угла между проекциями
    return np.degrees(np.arctan2(y, x))

def compute_phi_psi(atoms,chain,res_num):
    curr_atoms=atoms_in_residue(atoms,chain,res_num)
    if not curr_atoms:
        return (None, None)
    N_cur=atoms_xyz(curr_atoms, "N")
    CA_cur=atoms_xyz(curr_atoms, "CA")
    C_cur=atoms_xyz(curr_atoms, "C")
    if N_cur is None or CA_cur is None or C_cur is None:
        return (None, None)

    phi=None
    prev_atoms= atoms_in_residue(atoms,chain,str(int(res_num)-1))
    if prev_atoms:
        C_prev= atoms_xyz(prev_atoms, 'C')
        if C_prev is not None:
            phi=dihedral(C_prev, N_cur, CA_cur, C_cur)

    psi=None
    next_atoms= atoms_in_residue(atoms,chain,str(int(res_num)+1))
    if next_atoms:
        N_next= atoms_xyz(next_atoms, 'N')
        if N_next is not None:
            psi=dihedral(N_cur, CA_cur, C_cur, N_next)

    return (phi, psi)


#if __name__ == '__main__':
#    if len(sys.argv) < 4:
#        print("Использование: python mutate.py <pdb_file> <chain> <resseq>")
#        sys.exit(1)
#
#    pdb_file = sys.argv[1]
#    chain = sys.argv[2]
#    resseq = sys.argv[3]
#
#    atoms = read_pdb(pdb_file)
#    phi, psi = compute_phi_psi(atoms, chain, resseq)
#
#    print(f"Остаток {chain}{resseq}:")
#    if phi is not None:
#        print(f"  phi = {phi:.1f}°")
#    else:
#        print("  phi = не определён")
#    if psi is not None:
#        print(f"  psi = {psi:.1f}°")
#    else:
#        print("  psi = не определён")



#БЛОК 2 ЗАГРУЗКА БИБЛИОТЕК ПОИСК ПОДХОДЯЩЕГО РОТАМЕРА



def load_rotamer_library(lib_path):
    lib = {}
    with open(lib_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            res = parts[0]
            phi_bin = int(parts[1])
            psi_bin = int(parts[2])
            chi1 = float(parts[5])  
            chi2 = float(parts[7])   
            chi3 = float(parts[9])   
            chi4 = float(parts[11])
        

            key = (res, phi_bin, psi_bin)
            prob = float(parts[3])
            rotamer = {
                'prob': prob,
                'chi1': chi1, 'chi2': chi2, 'chi3': chi3, 'chi4': chi4
            }
            
            lib.setdefault(key, []).append(rotamer)
    return lib


def nearest_bin(angle): #округляю угол
    if angle is None:
        return None
    bin_angle=round(angle/10)*10
    
    candidate1 = bin_angle - 5
    candidate2 = bin_angle + 5
    if abs(angle - candidate1) <= abs(angle - candidate2):
        bin_angle = candidate1
    else:
        bin_angle = candidate2

    if bin_angle>175:
        bin_angle=175
    elif bin_angle< -175:
        bin_angle=-175

    return bin_angle



def best_rotamer(lib, target_res, phi, psi):
    b_phi = nearest_bin(phi)
    b_psi = nearest_bin(psi)
    target_res = target_res.upper()
    key = (target_res, b_phi, b_psi)
    if key in lib:
        rots = lib[key]
    else:
        candidates = [k for k in lib.keys() if k[0] == target_res]
        if not candidates:
            return None
        def dist(k):
            return (k[1] - b_phi) ** 2 + (k[2] - b_psi) ** 2
        best_key = min(candidates, key=dist)
        rots = lib[best_key]
    best = max(rots, key=lambda r: r['prob'])
    return best

    



#if __name__ == '__main__':
#    if len(sys.argv) < 4:
#        print("Использование: python mutate.py <pdb_file> <chain> <resseq>")
#        sys.exit(1)
#    pdb_file = sys.argv[1]
#    chain = sys.argv[2]
#    resseq = sys.argv[3]
#
#    atoms = read_pdb(pdb_file)
#    phi, psi = compute_phi_psi(atoms, chain, resseq)
#
#    print(f"Остаток {chain}{resseq}:")
#    if phi is not None:
#        print(f"  phi = {phi:.1f}°")
#    else:
#        print("  phi = не определён")
#    if psi is not None:
#        print(f"  psi = {psi:.1f}°")
#    else:
#        print("  psi = не определён")
#
#    lib_path = "ALL.bbdep.rotamers.lib"
#    lib = load_rotamer_library(lib_path)
#    print(f"Библиотека загружена. Всего ключей: {len(lib)}")
#
#    if phi is not None and psi is not None:
#        rot = best_rotamer(lib, "ARG", phi, psi)
#        if rot:
#            print(f"Лучший ротамер для ARG при phi={phi:.1f}, psi={psi:.1f}:")
#            print(f"  prob={rot['prob']:.3f}, chi1={rot['chi1']:.1f}, chi2={rot['chi2']:.1f}, chi3={rot['chi3']:.1f}, chi4={rot['chi4']:.1f}")
#        else:
#            print("Ротамер не найден")
#    else:
#        print("Углы не определены, проверка ротамеров пропущена")





#КООРДИНАТЫ НОВОГО АТОМА

def add_atom_by_internal(prev3, prev2, prev1, bond_length, bond_angle_deg, dihedral_deg):
    ang = np.radians(bond_angle_deg) #между прев 1 и новым     
    dihe = np.radians(dihedral_deg) 
    #локальная система координат по трем известным атомам
    z_axis=prev1-prev2
    z_axis=z_axis/np.linalg.norm(z_axis)#ось по направлению вдоль сведи между 1 и 2
    x_axis = prev2 - prev3
    x_axis = x_axis - np.dot(x_axis, z_axis) * z_axis
    if np.linalg.norm(x_axis) < 1e-12:
            x_axis = np.array([1.0, 0.0, 0.0])
            x_axis = x_axis - np.dot(x_axis, z_axis) * z_axis
    x_axis = x_axis / np.linalg.norm(x_axis)
    y_axis = np.cross(z_axis, x_axis)

    dx = bond_length*np.sin(ang)*np.cos(dihe)
    dy = bond_length*np.sin(ang)*np.sin(dihe)
    dz = bond_length*np.cos(ang)

    new_coord = prev1+dx*x_axis+dy*y_axis+dz*z_axis
    return new_coord
    

#cтандартная геометрия аминокислот
BOND_LENGTHS = {
    ('CA', 'CB'): 1.53,
    ('CB', 'CG'): 1.53,
    ('CG', 'CD'): 1.53,
    ('CD', 'NE'): 1.47,
    ('NE', 'CZ'): 1.33,
    ('CZ', 'NH1'): 1.33,
    ('CZ', 'NH2'): 1.33,
    ('CB', 'OG'): 1.43,
    ('OG', 'HG'): 0.96,
}

BOND_ANGLES = {
    ('N', 'CA', 'CB'): 110.5,
    ('CA', 'CB', 'CG'): 113.5,
    ('CB', 'CG', 'CD'): 113.5,
    ('CG', 'CD', 'NE'): 111.0,
    ('CD', 'NE', 'CZ'): 124.0,
    ('NE', 'CZ', 'NH1'): 120.0,
    ('NE', 'CZ', 'NH2'): 120.0,
    ('CA', 'CB', 'OG'): 110.5,
}



#для каждого остатка список шагов построения.

SIDECHAIN_TOPOLOGY = {
    'GLY': [],   
    'ALA': [],  
    'SER': [
        ('OG', ('N', 'CA', 'CB'), 1.43, 110.5, 0),  
    ],
    'CYS': [
        ('SG', ('N', 'CA', 'CB'), 1.81, 110.5, 0),   
    ],
    'VAL': [
        ('CG1', ('N', 'CA', 'CB'), 1.53, 110.5, 0),  
        ('CG2', ('CA', 'CB', 'CG1'), 1.53, 110.5, 0), 
    ],
    'LEU': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 113.5, 0),   
        ('CD1', ('CA', 'CB', 'CG'), 1.53, 110.5, 1),  
        ('CD2', ('CA', 'CB', 'CG'), 1.53, 110.5, 1),  
    ],
    'ILE': [
        ('CG1', ('N', 'CA', 'CB'), 1.53, 110.5, 0),  
        ('CG2', ('CA', 'CB', 'CG1'), 1.53, 110.5, 0), 
        ('CD1', ('CB', 'CG1', 'CG2'), 1.53, 110.5, 1), 
    ],
    'MET': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 113.5, 0),   
        ('SD', ('CA', 'CB', 'CG'), 1.81, 110.5, 1),  
        ('CE', ('CB', 'CG', 'SD'), 1.79, 100.0, 2), 
    ],
    'PHE': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 113.5, 0),   
        ('CD1', ('CA', 'CB', 'CG'), 1.51, 120.0, 1),  
        ('CD2', ('CA', 'CB', 'CG'), 1.51, 120.0, 1),  
    ],
    'TYR': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 113.5, 0),   
        ('CD1', ('CA', 'CB', 'CG'), 1.51, 120.0, 1),
        ('CD2', ('CA', 'CB', 'CG'), 1.51, 120.0, 1),
    ],
    'TRP': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 113.5, 0),
        ('CD1', ('CA', 'CB', 'CG'), 1.51, 120.0, 1),
        ('CD2', ('CA', 'CB', 'CG'), 1.51, 120.0, 1),
    ],
    'HIS': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 113.5, 0),
        ('ND1', ('CA', 'CB', 'CG'), 1.38, 120.0, 1),
        ('CD2', ('CA', 'CB', 'CG'), 1.38, 120.0, 1),
    ],
    'ASP': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 113.5, 0),
        ('OD1', ('CA', 'CB', 'CG'), 1.25, 118.0, 1),
        ('OD2', ('CA', 'CB', 'CG'), 1.25, 118.0, 1),
    ],
    'ASN': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 113.5, 0),
        ('OD1', ('CA', 'CB', 'CG'), 1.23, 120.0, 1),
        ('ND2', ('CA', 'CB', 'CG'), 1.33, 116.0, 1),
    ],
    'GLU': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 113.5, 0),
        ('CD', ('CA', 'CB', 'CG'), 1.53, 113.5, 1),
        ('OE1', ('CB', 'CG', 'CD'), 1.25, 118.0, 2),
        ('OE2', ('CB', 'CG', 'CD'), 1.25, 118.0, 2),
    ],
    'GLN': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 113.5, 0),
        ('CD', ('CA', 'CB', 'CG'), 1.53, 113.5, 1),
        ('OE1', ('CB', 'CG', 'CD'), 1.23, 120.0, 2),
        ('NE2', ('CB', 'CG', 'CD'), 1.33, 116.0, 2),
    ],
    'LYS': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 113.5, 0),
        ('CD', ('CA', 'CB', 'CG'), 1.53, 113.5, 1),
        ('CE', ('CB', 'CG', 'CD'), 1.53, 111.0, 2),
        ('NZ', ('CG', 'CD', 'CE'), 1.47, 110.0, 3),
    ],
    'ARG': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 113.5, 0),
        ('CD', ('CA', 'CB', 'CG'), 1.53, 113.5, 1),
        ('NE', ('CB', 'CG', 'CD'), 1.47, 111.0, 2),
        ('CZ', ('CG', 'CD', 'NE'), 1.33, 124.0, 3),
    ],
    'PRO': [
        ('CG', ('N', 'CA', 'CB'), 1.53, 104.0, 0),
        ('CD', ('CA', 'CB', 'CG'), 1.50, 105.0, 1), 
    ],
}



def build_sidechain(res_name, N, CA, C, chi_list): #новая боковая цепь
    atoms=[]
    if res_name.upper()!="GLY":
        CB = add_atom_by_internal(C, N, CA,
                                  bond_length=BOND_LENGTHS.get(('CA', 'CB'), 1.53),
                                  bond_angle_deg=BOND_ANGLES.get(('N', 'CA', 'CB'), 110.5),
                                  dihedral_deg=122.55)
        atoms.append(('CB', *CB)) #первый атом боковой цепи
    else:
        return atoms

    coords = {'N': N, 'CA': CA, 'C': C, 'CB': CB}
    top = SIDECHAIN_TOPOLOGY.get(res_name.upper(), [])
    for step in top:
        name, (p3_name, p2_name, p1_name), length, angle, chi_idx = step
        p3 = coords[p3_name]
        p2 = coords[p2_name]
        p1 = coords[p1_name]
        dihe = chi_list[chi_idx] if (chi_idx >= 0 and chi_idx < len(chi_list)) else 0.0
        new = add_atom_by_internal(p3, p2, p1, length, angle, dihe)
        atoms.append((name, *new))
        coords[name] = new
    return atoms


def mutate_residue(atoms, chain, res_num, new_res_name, new_sidechain):
    res_str = str(res_num)
    first_idx = None
    last_idx = None
    for i, line in enumerate(atoms):
        if len(line) < 26:
            continue
        if line[21:22] == chain and line[22:26].strip() == res_str:
            if first_idx is None:
                first_idx = i
            last_idx = i
    if first_idx is None:
        return atoms[:]

    before = atoms[:first_idx]
    after = atoms[last_idx+1:]

    keep_names = {'N', 'CA', 'C', 'O'}
    kept_atoms = []
    for line in atoms[first_idx:last_idx+1]:
        name = line[12:16].strip()
        if name in keep_names:
            new_line = line[:17] + f"{new_res_name:3s}" + line[20:]
            kept_atoms.append(new_line)

    new_atom_lines = []
    for name, x, y, z in new_sidechain:
        element = ''.join(c for c in name if c.isalpha())[0]
        line = (f"ATOM  {99999:5d} {name:4s} {new_res_name:3s} "
                f"{chain}{res_str:>4s}    {x:8.3f}{y:8.3f}{z:8.3f}"
                f"  1.00  0.00          {element:2s}  ")
        new_atom_lines.append(line)

    combined = before + kept_atoms + new_atom_lines + after

    for i, line in enumerate(combined):
        serial = i + 1
        combined[i] = f"{line[:6]}{serial:5d}{line[11:]}"

    return combined

def write_pdb(atom_lines, output_path):
    with open(output_path, 'w') as f:
        for line in atom_lines:
            f.write(line.rstrip()+'\n')










if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Мутагенез с ротамерной библиотекой Данбрэка')
    parser.add_argument('pdb_file', help='Входной PDB-файл')
    parser.add_argument('chain', help='Цепь (например, A)')
    parser.add_argument('res_num', help='Номер остатка (например, 42)')
    parser.add_argument('--to', '-t', required=True, help='Новая аминокислота (например, SER)')
    parser.add_argument('--out', '-o', required=True, help='Выходной PDB-файл')
    parser.add_argument('--lib', default='ALL.bbdep.rotamers.lib', help='Путь к библиотеке')
    args = parser.parse_args()

    atoms = read_pdb(args.pdb_file)
    phi, psi = compute_phi_psi(atoms, args.chain, args.res_num)
    print(f"Остаток {args.chain}{args.res_num}: phi={phi:.1f}°, psi={psi:.1f}°")

    if phi is None or psi is None:
        print("Не удалось вычислить оба угла. Мутация невозможна.")
        sys.exit(1)

    lib = load_rotamer_library(args.lib)


    target_res = args.to.upper()
    rot = best_rotamer(lib, target_res, phi, psi)
    if rot is None:
        print(f"Для {target_res} при данных углах ротамер не найден.")
        sys.exit(1)

    print(f"Ротамер: prob={rot['prob']:.1f}%, chi1={rot['chi1']:.1f}°, chi2={rot['chi2']:.1f}°, chi3={rot['chi3']:.1f}°, chi4={rot['chi4']:.1f}°")

    curr_atoms = atoms_in_residue(atoms, args.chain, args.res_num)
    N = atoms_xyz(curr_atoms, 'N')
    CA = atoms_xyz(curr_atoms, 'CA')
    C = atoms_xyz(curr_atoms, 'C')
    if N is None or CA is None or C is None:
        print("Не удалось получить координаты N, CA, C.")
        sys.exit(1)

    
    chi_angles = [rot['chi1'], rot['chi2'], rot['chi3'], rot['chi4']]
    new_sidechain = build_sidechain(target_res, N, CA, C, chi_angles)

    mutated_atoms = mutate_residue(atoms, args.chain, args.res_num, target_res, new_sidechain)

    write_pdb(mutated_atoms, args.out)
    print(f"Мутация в {args.out}")

















