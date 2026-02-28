#!/usr/bin/env python3

import sys
import os
import argparse


############################################################
# Utilities
############################################################

def strip_comments(line):
    return line.split(";")[0].strip()


def tokenize(line):
    return strip_comments(line).split()


############################################################
# Parse atom types used in input topology
############################################################

def get_atomtypes_from_itp(itp_file):

    atomtypes = set()
    section = None

    with open(itp_file) as f:
        for line in f:

            clean = strip_comments(line)

            if clean.startswith("["):
                section = clean.lower()
                continue

            if section == "[ atoms ]":
                fields = tokenize(line)
                if len(fields) >= 2:
                    atomtypes.add(fields[1])

    return atomtypes


############################################################
# Parse CHARMM force field atomtypes
############################################################

def parse_ff_atomtypes(ffdir):

    ff_atomtypes = {}

    path = os.path.join(ffdir, "ffnonbonded.itp")
    if not os.path.exists(path):
        raise FileNotFoundError(
            "ffnonbonded.itp not found in force field directory"
        )

    section = None

    with open(path) as f:
        for line in f:

            clean = strip_comments(line)

            if clean.startswith("["):
                section = clean.lower()
                continue

            if section == "[ atomtypes ]":

                if not clean:
                    continue

                fields = tokenize(line)

                if len(fields) >= 7:
                    atype = fields[0]
                    ff_atomtypes[atype] = line.rstrip()

    return ff_atomtypes


############################################################
# Write output
############################################################

def write_atomtypes(output_file, used_types, ff_atomtypes, strict=False):

    missing = []

    with open(output_file, "w") as out:

        out.write("[ atomtypes ]\n")
        out.write("; Extracted atomtypes from CHARMM force field\n\n")

        for atype in sorted(used_types):

            if atype in ff_atomtypes:
                out.write(ff_atomtypes[atype] + "\n")
            else:
                missing.append(atype)

    if missing:
        message = "Missing atomtypes in force field: " + ", ".join(missing)
        if strict:
            raise RuntimeError(message)
        else:
            print("WARNING:", message)


############################################################
# CLI
############################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Extract used CHARMM atomtypes into standalone .itp"
    )

    parser.add_argument("input_itp")
    parser.add_argument("ffdir")
    parser.add_argument("output_itp")
    parser.add_argument("--strict", action="store_true")

    args = parser.parse_args()

    used_types = get_atomtypes_from_itp(args.input_itp)
    ff_atomtypes = parse_ff_atomtypes(args.ffdir)

    write_atomtypes(
        args.output_itp,
        used_types,
        ff_atomtypes,
        strict=args.strict
    )