#!/usr/bin/env python3

import sys
import os
from collections import defaultdict

############################################
# Utilities
############################################

def strip_comments(line):
    return line.split(";")[0].strip()


def tokenize(line):
    return strip_comments(line).split()


############################################
# Force field parsing
############################################

def parse_ffnonbonded(ffdir):
    atomtypes = {}

    path = os.path.join(ffdir, "ffnonbonded.itp")
    if not os.path.exists(path):
        raise FileNotFoundError("ffnonbonded.itp not found")

    section = None
    with open(path) as f:
        for line in f:
            clean = strip_comments(line)
            if not clean:
                continue

            if clean.startswith("["):
                section = clean
                continue

            if section == "[ atomtypes ]":
                fields = tokenize(line)
                if len(fields) >= 7:
                    atomtypes[fields[0]] = {
                        "mass": fields[2],
                        "charge": fields[3],
                        "sigma": fields[5],
                        "epsilon": fields[6]
                    }

    return atomtypes


def parse_ffbonded(ffdir):

    bonds = {}
    angles = {}
    dihedrals = defaultdict(list)

    path = os.path.join(ffdir, "ffbonded.itp")
    if not os.path.exists(path):
        raise FileNotFoundError("ffbonded.itp not found")

    section = None
    with open(path) as f:
        for line in f:
            clean = strip_comments(line)
            if not clean:
                continue

            if clean.startswith("["):
                section = clean
                continue

            fields = tokenize(line)

            ################################
            # Bondtypes
            ################################
            if section == "[ bondtypes ]" and len(fields) >= 5:
                t1, t2 = fields[0], fields[1]
                funct = fields[2]
                r0 = fields[3]
                kb = fields[4]
                bonds[(t1, t2)] = (funct, r0, kb)

            ################################
            # Angletypes (CHARMM harmonic)
            ################################
            elif section == "[ angletypes ]" and len(fields) >= 6:
                t1, t2, t3 = fields[0:3]
                funct = fields[3]
                theta0 = fields[4]
                ktheta = fields[5]
                angles[(t1, t2, t3)] = (funct, theta0, ktheta)

            ################################
            # Dihedraltypes (CHARMM periodic)
            ################################
            elif section == "[ dihedraltypes ]" and len(fields) >= 8:
                t1, t2, t3, t4 = fields[0:4]
                funct = fields[4]
                phi0 = fields[5]
                kphi = fields[6]
                mult = fields[7]
                dihedrals[(t1, t2, t3, t4)].append(
                    (funct, phi0, kphi, mult)
                )

    return bonds, angles, dihedrals


############################################
# ITP annotation
############################################

def annotate_itp(itp_file, ffdir, output_file):

    atomtypes = parse_ffnonbonded(ffdir)
    bondtypes, angletypes, dihedraltypes = parse_ffbonded(ffdir)

    atom_index_to_type = {}
    section = None

    with open(itp_file) as fin, open(output_file, "w") as fout:

        for line in fin:

            clean = strip_comments(line)

            if clean.startswith("["):
                section = clean
                fout.write(line)
                continue

            if not clean:
                fout.write(line)
                continue

            fields = clean.split()

            ################################
            # ATOMS
            ################################
            if section == "[ atoms ]" and len(fields) >= 2:
                atomnr = fields[0]
                atype = fields[1]
                atom_index_to_type[atomnr] = atype
                fout.write(line)
                continue

            ################################
            # BONDS
            ################################
            if section == "[ bonds ]" and len(fields) >= 2:

                i, j = fields[0], fields[1]

                if i in atom_index_to_type and j in atom_index_to_type:
                    t1 = atom_index_to_type[i]
                    t2 = atom_index_to_type[j]

                    key = (t1, t2)
                    key_rev = (t2, t1)

                    if key in bondtypes:
                        funct, r0, kb = bondtypes[key]
                    elif key_rev in bondtypes:
                        funct, r0, kb = bondtypes[key_rev]
                    else:
                        fout.write(line)
                        continue

                    fout.write(
                        f"{i} {j} {funct} {r0} {kb}\n"
                    )
                    continue

            ################################
            # ANGLES
            ################################
            if section == "[ angles ]" and len(fields) >= 3:

                i, j, k = fields[0:3]

                if all(x in atom_index_to_type for x in (i, j, k)):
                    t1 = atom_index_to_type[i]
                    t2 = atom_index_to_type[j]
                    t3 = atom_index_to_type[k]

                    key = (t1, t2, t3)
                    key_rev = (t3, t2, t1)

                    if key in angletypes:
                        funct, theta0, ktheta = angletypes[key]
                    elif key_rev in angletypes:
                        funct, theta0, ktheta = angletypes[key_rev]
                    else:
                        fout.write(line)
                        continue

                    fout.write(
                        f"{i} {j} {k} {funct} {theta0} {ktheta}\n"
                    )
                    continue

            ################################
            # DIHEDRALS
            ################################
            if section == "[ dihedrals ]" and len(fields) >= 4:

                i, j, k, l = fields[0:4]

                if all(x in atom_index_to_type for x in (i, j, k, l)):

                    t1 = atom_index_to_type[i]
                    t2 = atom_index_to_type[j]
                    t3 = atom_index_to_type[k]
                    t4 = atom_index_to_type[l]

                    key = (t1, t2, t3, t4)
                    key_rev = (t4, t3, t2, t1)

                    terms = None
                    if key in dihedraltypes:
                        terms = dihedraltypes[key]
                    elif key_rev in dihedraltypes:
                        terms = dihedraltypes[key_rev]

                    if terms:
                        for funct, phi0, kphi, mult in terms:
                            fout.write(
                                f"{i} {j} {k} {l} {funct} {phi0} {kphi} {mult}\n"
                            )
                        continue

            fout.write(line)


############################################
# Main
############################################

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage:")
        print("python write_full_itp.py input.itp charmm36.ff output.itp")
        sys.exit(1)

    annotate_itp(sys.argv[1], sys.argv[2], sys.argv[3])
