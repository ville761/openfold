import argparse
import logging
import os
import sys
from pathlib import Path
import subprocess

from openfold.data.tools import hhsearch
from openfold.data.tools import mmseqs2

def _split_a3ms(output_dir):
    for fname in os.listdir(output_dir):
        if(not os.path.splitext(fname)[-1] == ".a3m"):
            continue

        fpath = os.path.join(output_dir, fname)

        with open(fpath, "r") as fp:
            a3ms = fp.read()

        # Split by the null byte, excluding the terminating null byte
        a3ms = a3ms.split('\x00')[:-1]

        for a3m in a3ms:
            name = a3m.split('\n', 1)[0][1:]
            prot_dir = os.path.join(output_dir, name)
            Path(prot_dir).mkdir(parents=True, exist_ok=True)
            with open(os.path.join(prot_dir, fname), "w") as fp:
                fp.write(a3m)

        os.remove(fpath)
        os.remove(fpath + ".dbtype")
        os.remove(fpath + ".index")


def main(args):
    with open(args.input_fasta, "r") as f:
        lines = [l.strip() for l in f.readlines()]

    names = lines[::2]
    seqs =  lines[1::2]

    if(args.fasta_chunk_size is None):
        chunk_size = len(seqs)
    else:
        chunk_size = args.fasta_chunk_size

    # Make the output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)


    s = 0
    while(s < len(seqs)):
        e = s + chunk_size
        chunk_fasta = [el for tup in zip(names[s:e], seqs[s:e]) for el in tup] 
        s = e
        
        #prot_dir = os.path.join(args.output_dir, chunk_fasta[0][1:].upper())
        prot_dir = os.path.join(args.output_dir, chunk_fasta[0][1:])
        if(os.path.exists(prot_dir)):
            # We've already computed this chunk
            continue

        mmseqs2.run_mmseqs(chunk_fasta[1],prot_dir)

        _split_a3ms(args.output_dir)


    hhsearch_pdb70_runner = hhsearch.HHSearch(
        binary_path=args.hhsearch_binary_path, databases=[args.pdb70]
    )


    for d in os.listdir(args.output_dir):
        dpath = os.path.join(args.output_dir, d)
        if(not os.path.isdir(dpath)):
            continue
        for fname in os.listdir(dpath):
            fpath = os.path.join(dpath, fname)
            if(not "uniref" in fname or 
                not os.path.splitext(fname)[-1] == ".a3m"):
                continue

            with open(fpath, "r") as fp:
                a3m = fp.read()

            hhsearch_result = hhsearch_pdb70_runner.query(a3m)
            pdb70_out_path = os.path.join(dpath, "pdb70_hits.hhr")
            with open(pdb70_out_path, "w") as f:
                f.write(hhsearch_result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_fasta", type=str,
        help="Path to input FASTA file. Can contain one or more sequences."
    )
    parser.add_argument(
        "uniref_db", type=str,
        help="Basename of uniref database"
    )
    parser.add_argument(
        "output_dir", type=str,
        help="Output directory"
    )
    parser.add_argument(
        "--hhsearch_binary_path", type=str, default=None,
        help="""Path to hhsearch binary (for template search). In future 
                versions, we'll also use mmseqs for this"""
    )
    parser.add_argument(
        "--pdb70", type=str, default=None,
        help="Basename of the pdb70 database"
    )
    parser.add_argument(
        "--env_db", type=str, default=None,
        help="Basename of environmental database"
    )
    parser.add_argument(
        "--fasta_chunk_size", type=int, default=None,
        help="""How many sequences should be processed at once. All sequences 
                processed at once by default."""
    )

    args = parser.parse_args()

    if(args.hhsearch_binary_path is not None and args.pdb70 is None):
        raise ValueError(
            "pdb70 must be specified along with hhsearch_binary_path"
        )

    main(args)
