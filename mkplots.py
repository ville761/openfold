#!/usr/bin/env python

import sys
import pickle
import os
import numpy as np
import argparse
import json
import matplotlib.pyplot as plt

def plot_ticks(Ls):
    Ln = sum(Ls)
    L_prev = 0
    for L_i in Ls[:-1]:
        L = L_prev + L_i
        L_prev += L_i
        plt.plot([0, Ln], [L, L], color="black")
        plt.plot([L, L], [0, Ln], color="black")
    ticks = np.cumsum([0]+Ls)
    ticks = (ticks[1:] + ticks[:-1])/2
    plt.yticks(ticks, alphabet_list[:len(ticks)])


def plot_plddts(plddts, Ls=None, dpi=100, fig=True):
    if fig:
        plt.figure(figsize=(8, 5), dpi=100)
    plt.title("Predicted lDDT per position")
    for n, plddt in enumerate(plddts):
        plt.plot(plddt, label=f"rank_{n+1}")
    if Ls is not None:
        L_prev = 0
        for L_i in Ls[:-1]:
            L = L_prev + L_i
            L_prev += L_i
            plt.plot([L, L], [0, 100], color="black")
    plt.legend()
    plt.ylim(0, 100)
    plt.ylabel("Predicted lDDT")
    plt.xlabel("Positions")
    return plt


def plot_paes(paes, Ls=None, dpi=100, fig=True):
    num_models = len(paes)
    if fig:
        plt.figure(figsize=(3*num_models, 2), dpi=dpi)
    for n, pae in enumerate(paes):
        plt.subplot(1, num_models, n+1)
        plt.title(f"rank_{n+1}")
        Ln = pae.shape[0]
        plt.imshow(pae, cmap="bwr", vmin=0, vmax=30, extent=(0, Ln, Ln, 0))
        if Ls is not None and len(Ls) > 1:
            plot_ticks(Ls)
        plt.colorbar()
    return plt


def plot_adjs(adjs, Ls=None, dpi=100, fig=True):
    num_models = len(adjs)
    if fig:
        plt.figure(figsize=(3*num_models, 2), dpi=dpi)
    for n, adj in enumerate(adjs):
        plt.subplot(1, num_models, n+1)
        plt.title(f"rank_{n+1}")
        Ln = adj.shape[0]
        plt.imshow(adj, cmap="binary", vmin=0, vmax=1, extent=(0, Ln, Ln, 0))
        if Ls is not None and len(Ls) > 1:
            plot_ticks(Ls)
        plt.colorbar()
    return plt


def plot_dists(dists, Ls=None, dpi=100, fig=True):
    num_models = len(dists)
    if fig:
        plt.figure(figsize=(3*num_models, 2), dpi=dpi)
    for n, dist in enumerate(dists):
        plt.subplot(1, num_models, n+1)
        plt.title(f"rank_{n+1}")
        Ln = dist.shape[0]
        plt.imshow(dist, extent=(0, Ln, Ln, 0))
        if Ls is not None and len(Ls) > 1:
            plot_ticks(Ls)
        plt.colorbar()
    return plt

def main(out,oudir,model_name,Lseq):

    plddt = out["plddt"]
    plddt_plot = plot_plddts([plddt], Ls=[Lseq], dpi=200)
    #plddt_png = args.output_dir.joinpath(tag + "_plddt.png")
    plddt_png = os.path.join(oudir,f"{model_name}_plddt.png")
    plddt_plot.savefig(str(plddt_png))
    plddt_plot.close()

    if "predicted_aligned_error" in out:
        pae_output_path = os.path.join(oudir, f'{model_name}_pae.json')
        pae = out["predicted_aligned_error"]
        max_pae = out['max_predicted_aligned_error']
        # Save predicted aligned error in the same format as the AF EMBL DB
        rounded_errors = np.round(pae.astype(np.float64), decimals=1)
        indices = np.indices((len(rounded_errors), len(rounded_errors))) + 1
        indices_1 = indices[0].flatten().tolist()
        indices_2 = indices[1].flatten().tolist()
        pae_data = json.dumps([{'residue1': indices_1,
                                'residue2': indices_2,
                                'distance': rounded_errors.flatten().tolist(),
                                'max_predicted_aligned_error': max_pae.item()
                              }],
                              indent=None,
                              separators=(',', ':'))
        with open(pae_output_path, 'w') as f:
            f.write(pae_data)

        paes_plot = plot_paes([pae], Ls=[Lseq], dpi=200)
        pae_png = os.path.join(oudir, f"{model_name}_PAE.png")
        paes_plot.savefig(str(pae_png))
        paes_plot.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "out_pickle", type=str,
    )
    parser.add_argument(
        "--model_name", type=str, required = True,
        help="""Model name""",
    )
    parser.add_argument(
        "--lseq", type=int, required = True,
        help="""sequence length""",
    )
    args = parser.parse_args()

    pf = args.out_pickle
    oudir = os.path.dirname(pf)
    if not oudir:
        oudir = "."

    try:
        with open(pf,'rb') as f:
            out = pickle.load(f)
    except IOError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        sys.exit(1)
    except:
        print("Unexpected error:", sys.exc_info()[0])
        sys.exit(1)

    main(out,oudir,args.model_name,args.lseq)
