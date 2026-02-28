#!/usr/bin/env python3

import sys
import os
from collections import defaultdict

###########################################################
# Utility functions
###########################################################

def strip_comments(line):
    return line.split(";")[0].strip()

def tokenize(line):
    return strip_comments(line).split()

def reverse_tuple(t):
    return tuple(reversed(t))

###########################################################
# Force field parsing
###########################################################

def parse_ffbonded(ffdir):

    bonds = []
    angles = []
    dihedrals = []
    impropers = []

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

            #################################################
            # BONDS
            #################################################
            if section == "[ bondtypes ]" and len(fields) >= 5:
                bonds.append({
                    "types": (fields[0], fields[1]),
                    "funct": fields[2],
                    "params": fields[3:5]
                })

            #################################################
            # ANGLES
            #################################################
            elif section == "[ angletypes ]" and len(fields) >= 6:
                angles.append({
                    "types": (fields[0], fields[1], fields[2]),
                    "funct": fields[3],
                    "params": fields[4:6]
                })

            #################################################
            # DIHEDRALS
            #################################################
            elif section == "[ dihedraltypes ]" and len(fields) >= 8:
                types = (fields[0], fields[1], fields[2], fields[3])
                funct = fields[4]
                params = fields[5:8]

                entry = {
                    "types": types,
                    "funct": funct,
                    "params": params
                }

                # CHARMM proper vs improper
                if funct == "2":  # harmonic improper
                    impropers.append(entry)
                else:
                    dihedrals.append(entry)

    return bonds, angles, dihedrals, impropers

###########################################################
# Wildcard matching with CHARMM specificity ranking
###########################################################

def match_types(query, candidate):

    wildcard_count = 0

    for q, c in zip(query, candidate):
        if c == "X":
            wildcard_count += 1
            continue
        if q != c:
            return None

    return wildcard_count


def resolve_parameter(query_types, database):

    best_match = None
    best_score = 999

    for entry in database:

        for candidate in (entry["types"], reverse_tuple(entry["types"])):

            score = match_types(query_types, candidate)

            if score is not None and score < best_score:
                best_match = entry
                best_score = score

    return best_match

###########################################################
# Main annotation logic
###########################################################

def annotate_itp(itp_file, ffdir, output_file):

    bonds_db, angles_db, dihedrals_db, impropers_db = parse_ffbonded(ffdir)

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

            #################################################
            # ATOMS
            #################################################
            if section == "[ atoms ]" and len(fields) >= 2:
                atom_index_to_type[fields[0]] = fields[1]
                fout.write(line)
                continue

            #################################################
            # BONDS
            #################################################
            if section == "[ bonds ]" and len(fields) >= 2:

                i, j = fields[0], fields[1]

                if i in atom_index_to_type and j in atom_index_to_type:

                    types = (
                        atom_index_to_type[i],
                        atom_index_to_type[j]
                    )

                    match = resolve_parameter(types, bonds_db)

                    if match:
                        funct = match["funct"]
                        r0, kb = match["params"]
                        fout.write(
                            f"{i} {j} {funct} {r0} {kb}\n"
                        )
                        continue

            #################################################
            # ANGLES
            #################################################
            if section == "[ angles ]" and len(fields) >= 3:

                i, j, k = fields[0:3]

                if all(x in atom_index_to_type for x in (i, j, k)):

                    types = (
                        atom_index_to_type[i],
                        atom_index_to_type[j],
                        atom_index_to_type[k]
                    )

                    match = resolve_parameter(types, angles_db)

                    if match:
                        funct = match["funct"]
                        theta0, ktheta = match["params"]
                        fout.write(
                            f"{i} {j} {k} {funct} {theta0} {ktheta}\n"
                        )
                        continue

            #################################################
            # DIHEDRALS
            #################################################
            if section == "[ dihedrals ]" and len(fields) >= 4:

                i, j, k, l = fields[0:4]

                if all(x in atom_index_to_type for x in (i, j, k, l)):

                    types = (
                        atom_index_to_type[i],
                        atom_index_to_type[j],
                        atom_index_to_type[k],
                        atom_index_to_type[l]
                    )

                    # Check improper first if funct == 2
                    if len(fields) >= 5 and fields[4] == "2":
                        db = impropers_db
                    else:
                        db = dihedrals_db

                    matches = []

                    for entry in db:
                        for candidate in (
                            entry["types"],
                            reverse_tuple(entry["types"])
                        ):
                            score = match_types(types, candidate)
                            if score is not None:
                                matches.append((score, entry))

                    if matches:
                        matches.sort(key=lambda x: x[0])
                        best_score = matches[0][0]

                        best_entries = [
                            e for s, e in matches if s == best_score
                        ]

                        for entry in best_entries:
                            funct = entry["funct"]
                            phi0, kphi, mult = entry["params"]
                            fout.write(
                                f"{i} {j} {k} {l} {funct} {phi0} {kphi} {mult}\n"
                            )
                        continue

            fout.write(line)


###########################################################
# CLI
###########################################################

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage:")
        print("python charmm_resolve.py input.itp charmm36.ff output.itp")
        sys.exit(1)

    annotate_itp(sys.argv[1], sys.argv[2], sys.argv[3])