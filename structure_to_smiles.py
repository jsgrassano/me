from typing import Optional, Iterable
import tempfile

# create an environment with:
# mamba env create -n smiles
# mamba install openbabel rdkit numpy
from rdkit import Chem

# needed in order to import pybel
import openbabel  # noqa
from openbabel import pybel

import numpy as np


def structure_to_raw_smiles(
    symbols: Iterable[str],
    coordinates: Iterable[float],
) -> Optional[str]:
    symbols = np.asarray(symbols).astype(str)
    coords = np.asarray(coordinates, dtype=np.float64)
    with tempfile.NamedTemporaryFile("w+") as f:
        f.write(f"{len(symbols)}\n")
        f.write("\n")
        for j, el in enumerate(symbols):
            f.write(f"{el} {coords[j][0]:8.3} {coords[j][1]:8.3} {coords[j][2]:8.3}\n")
        f.seek(0)
        obabel_molecule = next(pybel.readfile("xyz", f.name))
        raw_smiles = obabel_molecule.write(format="smi").split()[0].strip()
        # Not sure if raw_smiles can be None
    return raw_smiles


# RDKit returns None if it couldn't parse SMILES If the smiles can't be
# converted to an RDKit molecule the usual problem is hypervalent atoms. Since
# this happens for a lot of molecules, by default I perform only partial
# sanitization of the molecule and I don't enforce "traditional" valences.
#
# This partial sanitization procedure was found in the Rdkit mailing list, in
# [Rdkit-discuss] valence problem From: Adrian Jasiński <jasinski.adrian@gm...>
# - 2014-07-10 10:42:23 Attachments: Message as HTML link:
# https://sourceforge.net/p/rdkit/mailman/message/32589379/
_SANITIZATION_FLAGS = (
    Chem.SanitizeFlags.SANITIZE_FINDRADICALS
    | Chem.SanitizeFlags.SANITIZE_KEKULIZE
    | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
    | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
    | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
    | Chem.SanitizeFlags.SANITIZE_SYMMRINGS
)


def canonicalize_smiles(smiles: Optional[str]) -> Optional[str]:
    if smiles is None:
        return smiles
    assert smiles is not None  # mypy
    canonical_smiles = None
    try:
        molecule = Chem.MolFromSmiles(smiles, sanitize=True)
        if molecule is not None:
            canonical_smiles = Chem.MolToSmiles(molecule)
        else:
            molecule = Chem.MolFromSmiles(smiles, sanitize=False)
            # Only sanitize if the conversion worked
            if molecule is not None:
                molecule.UpdatePropertyCache(strict=False)
                Chem.SanitizeMol(molecule, _SANITIZATION_FLAGS, catchErrors=True)
                # Only continue if everything before worked
                if molecule:
                    canonical_smiles = Chem.MolToSmiles(molecule)
    except Exception as ex:
        print(ex)
        pass
    finally:
        return canonical_smiles


def structure_to_smiles(
    symbols: Iterable[str],
    coordinates: Iterable[float],
) -> Optional[str]:
    # if the function returns None, it could not create a smiles :'(
    return canonicalize_smiles(structure_to_raw_smiles(symbols, coordinates))









#if __name__ == "__main__":
#    from pathlib import Path
    # Example loading xyz files
#    for xyz in sorted(Path("/path/to/xyz/dir").resolve().iterdir()):
#        if xyz.suffix != ".xyz":
#            continue
#        symbols = []
#        coords = []
#        with open(xyz, "r") as f:
#            lines = iter(f.readlines())
#            next(lines)
#            next(lines)
#            for line in lines:
#                el, x, y, z = line.split()
#                symbols.append(el)
#                coords.append([float(x), float(y), float(z)])
#        smiles = structure_to_smiles(symbols, coords)
#        print(smiles)
