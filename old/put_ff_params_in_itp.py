#!/usr/bin/env python3

import sys
import os
import re
from collections import defaultdict

############################################
# Utility parsing helpers
############################################

def strip_comments(line):
    if ";" in line:
        return line.split(";")[0].strip()
    return line.strip()


def tokenize(line):
    return strip_comments(line).split()


############################################
# Force field parsing
############################################

def parse_ffnonbonded(ffdir):
    atomtypes = {}

    path = os.path.join(ffdir, "ffnonbonded.itp")
    if not os.path.exists(path):
        raise FileNotFoundError("ffnonbonded.itp not found in ff directory")

    section = None
    with open(path) as f:
        for line in f:
            line_clean = strip_comments(line)
            if not line_clean:
                continue

            if line_clean.startswith("["):
                section = line_clean.strip()
                continue

            if section == "[ atomtypes ]":
                fields = tokenize(line)
                if len(fields) >= 6:
                    atype = fields[0]
                    mass = fields[2]
                    charge = fields[3]
                    ptype = fields[4]
                    sigma = fields[5]
                    epsilon = fields[6] if len(fields) > 6 else ""
                    atomtypes[atype] = {
                        "mass": mass,
                        "charge": charge,
                        "ptype": ptype,
                        "sigma": sigma,
                        "epsilon": epsilon
                    }

    return atomtypes


def parse_ffbonded(ffdir):
    bonds = {}
    angles = {}
    dihedrals = defaultdict(list)

    path = os.path.join(ffdir, "ffbonded.itp")
    if not os.path.exists(path):
        raise FileNotFoundError("ffbonded.itp not found in ff directory")

    section = None
    with open(path) as f:
        for line in f:
            line_clean = strip_comments(line)
            if not line_clean:
                continue

            if line_clean.startswith("["):
                section = line_clean.strip()
                continue

            fields = tokenize(line)

            if section == "[ bondtypes ]" and len(fields) >= 5:
                key = tuple(fields[0:2])
                bonds[key] = fields[2:]

            elif section == "[ angletypes ]" and len(fields) >= 6:
                key = tuple(fields[0:3])
                angles[key] = fields[3:]

            elif section == "[ dihedraltypes ]" and len(fields) >= 8:
                key = tuple(fields[0:4])
                dihedrals[key].append(fields[4:])

    return bonds, angles, dihedrals


############################################
# ITP processing
############################################

def annotate_itp(itp_file, ffdir, output_file):

    atomtypes = parse_ffnonbonded(ffdir)
    bondtypes, angletypes, dihedraltypes = parse_ffbonded(ffdir)

    atom_index_to_type = {}

    section = None

    with open(itp_file) as fin, open(output_file, "w") as fout:

        for line in fin:
            stripped = strip_comments(line)

            if stripped.startswith("["):
                section = stripped
                fout.write(line)
                continue

            if not stripped:
                fout.write(line)
                continue

            fields = stripped.split()

            ####################################
            # ATOMS
            ####################################
            if section == "[ atoms ]":
                if len(fields) >= 2:
                    atomnr = fields[0]
                    atype = fields[1]
                    atom_index_to_type[atomnr] = atype

                    if atype in atomtypes:
                        p = atomtypes[atype]
                        annotation = (
                            f" ; sigma={p['sigma']} epsilon={p['epsilon']} "
                            f"mass={p['mass']} charge={p['charge']}"
                        )
                        fout.write(line.rstrip() + annotation + "\n")
                        continue

            ####################################
            # BONDS
            ####################################
            elif section == "[ bonds ]":
                if len(fields) >= 2:
                    i, j = fields[0], fields[1]
                    if i in atom_index_to_type and j in atom_index_to_type:
                        t1 = atom_index_to_type[i]
                        t2 = atom_index_to_type[j]
                        key = (t1, t2)
                        key_rev = (t2, t1)

                        if key in bondtypes:
                            params = bondtypes[key]
                        elif key_rev in bondtypes:
                            params = bondtypes[key_rev]
                        else:
                            params = None

                        if params:
                            annotation = " ; " + " ".join(params)
                            fout.write(line.rstrip() + annotation + "\n")
                            continue

            ####################################
            # ANGLES
            ####################################
            elif section == "[ angles ]":
                if len(fields) >= 3:
                    i, j, k = fields[0:3]
                    if all(x in atom_index_to_type for x in (i, j, k)):
                        t1 = atom_index_to_type[i]
                        t2 = atom_index_to_type[j]
                        t3 = atom_index_to_type[k]
                        key = (t1, t2, t3)
                        key_rev = (t3, t2, t1)

                        if key in angletypes:
                            params = angletypes[key]
                        elif key_rev in angletypes:
                            params = angletypes[key_rev]
                        else:
                            params = None

                        if params:
                            annotation = " ; " + " ".join(params)
                            fout.write(line.rstrip() + annotation + "\n")
                            continue

            ####################################
            # DIHEDRALS
            ####################################
            elif section == "[ dihedrals ]":
                if len(fields) >= 4:
                    i, j, k, l = fields[0:4]
                    if all(x in atom_index_to_type for x in (i, j, k, l)):
                        t1 = atom_index_to_type[i]
                        t2 = atom_index_to_type[j]
                        t3 = atom_index_to_type[k]
                        t4 = atom_index_to_type[l]

                        key = (t1, t2, t3, t4)
                        key_rev = (t4, t3, t2, t1)

                        params_list = None
                        if key in dihedraltypes:
                            params_list = dihedraltypes[key]
                        elif key_rev in dihedraltypes:
                            params_list = dihedraltypes[key_rev]

                        if params_list:
                            annotation = " ; "
                            for p in params_list:
                                annotation += "[" + " ".join(p) + "] "
                            fout.write(line.rstrip() + annotation + "\n")
                            continue

            fout.write(line)


############################################
# Main
############################################

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage:")
        print("python annotate_itp.py input.itp charmm36.ff output.itp")
        sys.exit(1)

    itp_file = sys.argv[1]
    ffdir = sys.argv[2]
    output_file = sys.argv[3]

    annotate_itp(itp_file, ffdir, output_file)

    print(f"Annotated topology written to {output_file}")