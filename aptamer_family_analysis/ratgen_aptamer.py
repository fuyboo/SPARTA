import argparse
import sys
import os
import logging
import random
import pandas as pd
import numpy as np
import torch
import torch.optim as optim
from pathlib import Path
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader
from collections import Counter
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import pysam


RAPGEN_MODELS_PATH = "/home/disk/fuyongbo/lgy/aptamer_pred/aptamer_analysis/Ratgen"
sys.path.append(RAPGEN_MODELS_PATH)
import raptgen_models

# Function to prepare training data
def prepare_data(ptk7_sample_df, negatibe_sample_df, output_path="raptgen_data/gen_train_data.fasta", sample_len=46):
    obj_sample_seq = [seq for seq in ptk7_sample_df['sequence'] if len(seq) == sample_len and 'N' not in seq]
    background_sample_seq = [seq for seq in negatibe_sample_df['sequence'] if len(seq) == sample_len and 'N' not in seq]
    
    sample_total_len = min(len(obj_sample_seq), len(background_sample_seq))
    random.seed(42)
    sampled_obj_sample_seq = random.sample(obj_sample_seq, sample_total_len)
    sampled_background_sample_seq = random.sample(background_sample_seq, sample_total_len)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f_io:
        for count, obj_seq in enumerate(sampled_obj_sample_seq):
            print(f'>pos_gen_train_seq_{count}', file=f_io)
            print(obj_seq, file=f_io)
        for count, obj_seq in enumerate(sampled_background_sample_seq):
            print(f'>nega_gen_train_seq_{count}', file=f_io)
            print(obj_seq, file=f_io)
    
    print(f"Data saved to {output_path}")

# Function to train the model
def train_model(data_path, save_dir, epochs=1000, batch_size=512, test_size=0.2, target_len=46, device=None):
    os.makedirs(save_dir, exist_ok=True)
    save_dir = Path(save_dir)
    
    # Read and process data
    c = raptgen_models.read_fasta(data_path)
    sequences = list(filter(lambda seq_count: seq_count[1] >= 1, Counter(c).most_common()))
    sequences = list(filter(lambda sequence: len(sequence) == target_len, sequences))
    seq, _ = zip(*sequences)
    
    train_test = np.array(list(map(raptgen_models.one_hot_index, seq)))
    train_data, test_data = train_test_split(train_test, test_size=test_size, shuffle=True, random_state=42)
    
    train_data = torch.from_numpy(train_data).long()
    test_data = torch.from_numpy(test_data).long()
    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_data, batch_size=batch_size, shuffle=False)

    # Model setup
    model = raptgen_models.CNN_PHMM_VAE(target_len, embed_size=2)
    model_str = "cnn_phmm_vae.mdl"
    optimizer = optim.Adam(model.parameters())
    model = model.to(device)

    train_kwargs = {
        "epochs": epochs,
        "train_loader": train_loader,
        "test_loader": test_loader,
        "save_dir": save_dir,
        "model": model,
        "model_str": model_str,
        "optimizer": optimizer,
        "device": device,
        "threshold": 50,
        "logger": logging.getLogger(__name__)
    }

    # Train the model
    logging.info(f"Training model: {model_str}")
    raptgen_models.train(**train_kwargs)

# Function to postprocess the model and visualize results
def postprocess_model(model_path, data_path, save_dir, eval_max=256, device=None):
    # Load the trained model
    model = raptgen_models.CNN_PHMM_VAE(46, embed_size=2)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model = model.to(device)

    # Encode sequences
    c = raptgen_models.read_fasta(data_path)
    sequences = list(filter(lambda seq_count: seq_count[1] >= 1, Counter(c).most_common()))
    sequences = list(filter(lambda sequence: len(sequence) == 46, sequences))
    seq, _ = zip(*sequences)

    mus = []
    with torch.no_grad():
        model.eval()
        for sequence in seq:
            recon, mu, logvar = model(torch.Tensor([raptgen_models.one_hot_index(sequence)]).long().to(device))
            mus += [mu]
    embed_x = torch.cat(mus)
    
    dims = embed_x.shape[1]
    with open(save_dir / "embed_seq.csv", "w") as f:
        f.write("index,seq," + ",".join([f"dim{dim+1}" for dim in range(dims)]) + "\n")
        for i, (seq, X) in enumerate(zip(sequences, embed_x)):
            f.write(f"{i},{seq[0]}," + ",".join(list(map(lambda x: f"{x}", X))) + "\n")

    # Clustering
    label_data_dict = {'name': [], 'seq': [], 'label': []}
    with pysam.FastxFile(data_path) as fq_io:
        for read in fq_io:
            label_data_dict['name'].append(read.name)
            label_data_dict['seq'].append(read.sequence)
            if read.name.split('_')[0] == 'pos':
                label_data_dict['label'].append(1)
            else:
                label_data_dict['label'].append(0)
    label_data_df = pd.DataFrame(label_data_dict).set_index('seq')
    embed_df = pd.read_csv(save_dir / "embed_seq.csv", index_col=0).set_index('seq')
    
    merged_df = pd.merge(embed_df, label_data_df, left_index=True, right_index=True, how='inner')
    
    plt.figure(figsize=(4, 4))
    plt.scatter(merged_df[merged_df['label'] == 1]['dim1'], merged_df[merged_df['label'] == 1]['dim2'], s=4, label='PTK7_apt')
    plt.scatter(merged_df[merged_df['label'] != 1]['dim1'], merged_df[merged_df['label'] != 1]['dim2'], s=4, label='Background')
    plt.legend()
    plt.title('VAE Latent Space')
    plt.savefig(save_dir / 'read_data_mix_label.jpg', bbox_inches='tight')
    merged_df.to_csv(save_dir / 'read_data_mix_label.csv')

# Main function to run the whole pipeline
def run_aptamer_training(ptk7_sample_path, negatibe_sample_path, model_save_dir='raptgen_save', data_save_dir='raptgen_data', device='cuda:3'):
    # Load the input data
    ptk7_sample_df = pd.read_csv(ptk7_sample_path, index_col=0)
    negatibe_sample_df = pd.read_csv(negatibe_sample_path, index_col=0)
    
    # Step 1: Prepare data
    gen_train_data_filename = Path(data_save_dir) / 'gen_train_data.fasta'
    prepare_data(ptk7_sample_df, negatibe_sample_df, output_path=gen_train_data_filename)
    
    # Step 2: Train model
    train_model(gen_train_data_filename, save_dir=model_save_dir, device=torch.device(device))
    
    # Step 3: Post-process model and visualize results
    model_path = Path(model_save_dir) / 'cnn_phmm_vae.mdl'
    postprocess_model(model_path, gen_train_data_filename, save_dir=model_save_dir, device=torch.device(device))

# Command-line interface
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aptamer training pipeline")
    parser.add_argument('-ptk7_sample_path', required=True, help="Path to PTK7 sample file")
    parser.add_argument('-negatibe_sample_path', required=True, help="Path to negative sample file")
    parser.add_argument('-model_save_dir', default='raptgen_save', help="Directory to save the trained model")
    parser.add_argument('-data_save_dir', default='raptgen_data', help="Directory to save generated data")
    parser.add_argument('-device', default='cuda:3', help="Device to use for training (e.g., 'cuda:0' or 'cpu')")

    args = parser.parse_args()
    run_aptamer_training(
        ptk7_sample_path=args.ptk7_sample_path,
        negatibe_sample_path=args.negatibe_sample_path,
        model_save_dir=args.model_save_dir,
        data_save_dir=args.data_save_dir,
        device=args.device
    )
