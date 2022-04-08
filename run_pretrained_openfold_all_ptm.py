# Copyright 2021 AlQuraishi Laboratory
# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# VU: This version produces all pTM models by default and highlights the best
# model.

import argparse
import logging
import numpy as np
import os

import pickle
import random
import sys
import time
import torch
import json

from openfold.config import model_config
from openfold.data import templates, feature_pipeline, data_pipeline
from openfold.model.model import AlphaFold
from openfold.np import residue_constants, protein
import openfold.np.relax.relax as relax
from openfold.utils.import_weights import (
    import_jax_weights_,
)
from openfold.utils.tensor_utils import (
    tensor_tree_map,
)

from scripts.utils import add_data_args

model_names = ['model_1_ptm', 'model_2_ptm', 'model_3_ptm', 'model_4_ptm',
               'model_5_ptm']

logging.basicConfig(filename='run_pretrained_openfold.log', filemode='w',
                    level=logging.DEBUG)


def main(args):
    config = model_config(model_names[0])

    template_featurizer = templates.TemplateHitFeaturizer(
        mmcif_dir=args.template_mmcif_dir,
        max_template_date=args.max_template_date,
        max_hits=config.data.predict.max_templates,
        kalign_binary_path=args.kalign_binary_path,
        release_dates_path=args.release_dates_path,
        obsolete_pdbs_path=args.obsolete_pdbs_path
    )

    use_small_bfd = (args.bfd_database_path is None)

    data_processor = data_pipeline.DataPipeline(
        template_featurizer=template_featurizer,
    )

    output_dir_base = args.output_dir
    if not os.path.exists(output_dir_base):
        os.makedirs(output_dir_base)
    random_seed = args.data_random_seed
    if random_seed is None:
        random_seed = random.randrange(sys.maxsize)
    if(args.use_precomputed_alignments is None):
        alignment_dir = os.path.join(output_dir_base, "alignments")
    else:
        alignment_dir = args.use_precomputed_alignments
    feature_processor = feature_pipeline.FeaturePipeline(config.data)

    # Gather input sequences
    with open(args.fasta_path, "r") as fp:
        lines = [line.strip() for line in fp.readlines()]

    tags, seqs = lines[::2], lines[1::2]
    tags = [tag[1:] for tag in tags]

    for tag, seq in zip(tags, seqs):
        fasta_path = os.path.join(args.output_dir, "tmp.fasta")
        with open(fasta_path, "w") as fp:
            fp.write(f">{tag}\n{seq}")

        logging.info("Generating features...")
        local_alignment_dir = os.path.join(alignment_dir, tag)
        if(args.use_precomputed_alignments is None):
            if not os.path.exists(local_alignment_dir):
                os.makedirs(local_alignment_dir)

            alignment_runner = data_pipeline.AlignmentRunner(
                jackhmmer_binary_path=args.jackhmmer_binary_path,
                hhblits_binary_path=args.hhblits_binary_path,
                hhsearch_binary_path=args.hhsearch_binary_path,
                uniref90_database_path=args.uniref90_database_path,
                mgnify_database_path=args.mgnify_database_path,
                bfd_database_path=args.bfd_database_path,
                uniclust30_database_path=args.uniclust30_database_path,
                pdb70_database_path=args.pdb70_database_path,
                use_small_bfd=use_small_bfd,
                no_cpus=args.cpus,
            )
            alignment_runner.run(
                fasta_path, local_alignment_dir
            )

        feature_dict = data_processor.process_fasta(
            fasta_path=fasta_path, alignment_dir=local_alignment_dir
        )

        # Remove temporary FASTA file
        os.remove(fasta_path)

        plddts = {}
        pae_outputs = {}
        unrelaxed_proteins = {}

        for model_name in model_names:
            processed_feature_dict = feature_processor.process_features(
                feature_dict, mode='predict',
            )
            config = model_config(model_name)
            model = AlphaFold(config)
            model = model.eval()
            param_path = os.path.join(
                "openfold", "resources", "params",
                "params_" + model_name + ".npz"
            )
            import_jax_weights_(model, param_path, version=model_name)
            model = model.to(args.model_device)
            logging.info("Executing model "+model_name+"...")
            batch = processed_feature_dict
            with torch.no_grad():
                batch = {
                    k: torch.as_tensor(v, device=args.model_device) for k, v in
                    batch.items()
                }

                t = time.perf_counter()
                out = model(batch)
                logging.info(f"Inference time: {time.perf_counter() - t}")

            # Toss out the recycling dimensions --- we don't need them anymore
            batch = tensor_tree_map(lambda x: np.array(x[..., -1].cpu()),
                                    batch)
            out = tensor_tree_map(lambda x: np.array(x.cpu()), out)

            plddt = out["plddt"]

            if 'predicted_aligned_error' in out:
                pae_outputs[model_name] = (
                                out['predicted_aligned_error'],
                                out['max_predicted_aligned_error']
                            )
            # not sure if this should be done for the pTM models. See the
            # OpenFold jupyter notebook.
            plddts[model_name] = plddt

            plddt_b_factors = np.repeat(
                plddt[..., None], residue_constants.atom_type_num, axis=-1
            )

            unrelaxed_protein = protein.from_prediction(
                features=batch,
                result=out,
                b_factors=plddt_b_factors
            )
            unrelaxed_proteins[model_name] = unrelaxed_protein

            # Save the model outputs.
            result_output_path = os.path.join(args.output_dir,
                                              f'result_{tag}_{model_name}.pkl')
            with open(result_output_path, 'wb') as f:
                pickle.dump(out, f, protocol=4)

            # Delete unused outputs to save memory.
            del model
            del processed_feature_dict
            del out

            # Save the unrelaxed PDB.
            unrelaxed_output_path = os.path.join(
                args.output_dir, f'{tag}_{model_name}_unrelaxed.pdb'
            )
            with open(unrelaxed_output_path, 'w') as f:
                f.write(protein.to_pdb(unrelaxed_protein))

        amber_relaxer = relax.AmberRelaxation(
                use_gpu=(args.model_device != "cpu"), **config.relax,)

        # Find the best model according to the mean pLDDT.
        best_model_name = max(plddts.keys(), key=lambda x: plddts[x].mean())

        # Relax the prediction.
        t = time.perf_counter()
        visible_devices = os.getenv("CUDA_VISIBLE_DEVICES")
        if("cuda" in args.model_device):
            device_no = args.model_device.split(":")[-1]
            os.environ["CUDA_VISIBLE_DEVICES"] = device_no
        relaxed_pdb_str, _, _ = amber_relaxer.process(
            prot=unrelaxed_proteins[best_model_name])
        if visible_devices:
            os.environ["CUDA_VISIBLE_DEVICES"] = visible_devices
        logging.info(f"Relaxation time: {time.perf_counter() - t}")

        # Save the relaxed PDB.
        relaxed_output_path = os.path.join(
                args.output_dir, f'{tag}_selected_prediction.pdb')
        with open(relaxed_output_path, 'w') as f:
            f.write(relaxed_pdb_str)

        # Save pLDDT and predicted aligned error (if it exists)
        pae_output_path = os.path.join(
                args.output_dir, 'predicted_aligned_error.json')
        if pae_outputs:
            pae, max_pae = pae_outputs[best_model_name]
            # Save predicted aligned error in the same format as the AF EMBL DB
            rounded_errors = np.round(pae.astype(np.float64), decimals=1)
            indices = np.indices(
                    (len(rounded_errors), len(rounded_errors))
                    ) + 1
            indices_1 = indices[0].flatten().tolist()
            indices_2 = indices[1].flatten().tolist()
            pae_data = json.dumps([{
                            'residue1': indices_1,
                            'residue2': indices_2,
                            'distance': rounded_errors.flatten().tolist(),
                            'max_predicted_aligned_error': max_pae.item()
                            }],
                            indent=None,
                            separators=(',', ':'))
            with open(pae_output_path, 'w') as f:
                f.write(pae_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fasta_path", type=str,
    )
    parser.add_argument(
        "template_mmcif_dir", type=str,
    )
    parser.add_argument(
        "--use_precomputed_alignments", type=str, default=None,
        help="""Path to alignment directory. If provided, alignment computation
                is skipped and database path arguments are ignored."""
    )
    parser.add_argument(
        "--output_dir", type=str, default=os.getcwd(),
        help="""Name of the directory in which to output the prediction""",
    )
    parser.add_argument(
        "--model_device", type=str, default="cpu",
        help="""Name of the device on which to run the model. Any valid torch
             device name is accepted (e.g. "cpu", "cuda:0")"""
    )
    parser.add_argument(
        "--cpus", type=int, default=4,
        help="""Number of CPUs with which to run alignment tools"""
    )
    parser.add_argument(
        '--preset', type=str, default='full_dbs',
        choices=('reduced_dbs', 'full_dbs')
    )
    parser.add_argument(
        '--data_random_seed', type=str, default=None
    )
    add_data_args(parser)
    args = parser.parse_args()

    if(args.model_device == "cpu" and torch.cuda.is_available()):
        logging.warning(
            """The model is being run on CPU. Consider specifying
            --model_device for better performance"""
        )

    main(args)
