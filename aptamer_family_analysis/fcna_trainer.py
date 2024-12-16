import os
import sys
import time
import argparse
import math
import numpy as np
import os.path as osp
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
import torch.nn.functional as F
import pandas as pd
import random
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import sklearn.metrics as metrics

# Define default parameters
weight_decay = 0.0005
learning_rate = 0.001
batch_size = 256
seq_dict = {'A': [1, 0, 0, 0], 'G': [0, 0, 1, 0],
            'C': [0, 1, 0, 0], 'T': [0, 0, 0, 1],
            'a': [1, 0, 0, 0], 'g': [0, 0, 1, 0],
            'c': [0, 1, 0, 0], 't': [0, 0, 0, 1],
            'N': [0, 0, 0, 0], 'n': [0, 0, 0, 0]}
seq_len = 48

# Function to upsample
def upsample(x, out_size):
    return F.interpolate(x, size=out_size, mode='linear', align_corners=False)

# Batch normalization, ReLU, and Convolution layers
def bn_relu_conv(in_, out_, kernel_size=3, stride=1, bias=False):
    padding = kernel_size // 2
    return nn.Sequential(nn.BatchNorm1d(in_),
                         nn.ReLU(inplace=True),
                         nn.Conv1d(in_, out_, kernel_size=kernel_size, stride=stride, padding=padding, bias=bias))

# FCNA model class
class FCNA(nn.Module):
    """FCNA for motif mining"""
    def __init__(self, motiflen=13):
        super(FCNA, self).__init__()
        self.conv1 = nn.Conv1d(in_channels=4, out_channels=64, kernel_size=motiflen)
        self.pool1 = nn.MaxPool1d(kernel_size=4, stride=4)
        self.conv2 = nn.Conv1d(in_channels=64, out_channels=64, kernel_size=5)
        self.pool2 = nn.MaxPool1d(kernel_size=4, stride=4)
        self.conv3 = nn.Conv1d(in_channels=64, out_channels=64, kernel_size=3)
        self.pool3 = nn.MaxPool1d(kernel_size=2, stride=2)
        self.aap = nn.AdaptiveAvgPool1d(1)
        # Decode layers
        self.blend4 = bn_relu_conv(64, 64, kernel_size=3)
        self.blend3 = bn_relu_conv(64, 64, kernel_size=3)
        self.blend2 = bn_relu_conv(64, 4, kernel_size=3)
        self.blend1 = bn_relu_conv(4, 1, kernel_size=3)
        self.linear1 = nn.Linear(seq_len, 2)
        # Activation, dropout
        self.sigmoid = nn.Sigmoid()
        self.relu = nn.ReLU(inplace=True)
        self.dropout = nn.Dropout(p=0.2)
        self._init_weights()

    def _init_weights(self):
        for layer in self.modules():
            if isinstance(layer, (nn.Conv1d, nn.Linear)):
                nn.init.xavier_uniform_(layer.weight)
                if layer.bias is not None:
                    nn.init.constant_(layer.bias, 0)
            elif isinstance(layer, nn.BatchNorm1d):
                nn.init.constant_(layer.weight, 1)
                nn.init.constant_(layer.bias, 0)

    def forward(self, data):
        b, _, _ = data.size()
        # Encoding process
        skip1 = data  
        out1 = self.conv1(data)
        out1 = self.relu(out1)
        out1 = self.pool1(out1)
        skip2 = out1
        out1 = self.conv2(out1)
        out1 = self.relu(out1)
        out1 = self.pool2(out1)
        out1 = self.dropout(out1)
        skip3 = out1
        up3 = self.aap(out1)
        up2 = upsample(up3, skip2.size()[-1])
        up2 = up2 + skip2
        up2 = self.blend2(up2)
        up1 = upsample(up2, skip1.size()[-1])
        up1 = up1 + skip1
        up1 = self.blend1(up1)
        out_dense = self.sigmoid(up1)
        out_dense = out_dense.view(b, -1)
        out_res = self.linear1(out_dense)
        return out_res

# Plotting confusion matrix
def plot_confusion_matrix(y_true, y_pred, classes, normalize=False, title=None, cmap=plt.cm.Blues, ax=None):
    plt.rcParams["figure.figsize"] = [6, 4]
    cm = metrics.confusion_matrix(y_true, y_pred, labels=classes)
    if normalize:
        cm = cm.astype("float") / cm.sum(axis=1)[:, np.newaxis]

    if ax is None:
        (fig, ax) = plt.subplots()

    im = ax.imshow(cm, interpolation="nearest", cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           xticklabels=classes,
           yticklabels=classes,
           title=title,
           ylabel="True label",
           xlabel="Predicted label")
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(12)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    fmt = ".2f" if normalize else "d"
    thresh = cm.max() / 2.0
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt), ha="center", va="center", color="white" if cm[i, j] > thresh else "black")

    return fig, ax

# Main function to run the training and evaluation
import pandas as pd
import torch
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
import random
import os

def main(args):
    external_data = args.external_data 
    train_data = args.train_data
    output_path = args.output_path
    os.makedirs(output_path, exist_ok=True)

    model = FCNA(motiflen=13)

    # Load data
    label_df = pd.read_csv(train_data, index_col=0)
    train_data_list, train_label_list, val_data_list, val_label_list = [], [], [], []

    random.seed(2021)
    val_list = random.sample(range(label_df.shape[0]), int(label_df.shape[0] * 0.2))
    train_list = list(set(range(label_df.shape[0])) - set(val_list))
    val_df = label_df.iloc[val_list]
    train_df = label_df.iloc[train_list]

    for idx, data in train_df.iterrows():
        need_line_str = str(data['sequence']).replace('\n', '')
        if len(need_line_str) <= seq_len:
            need_line_str = need_line_str.rjust(seq_len, 'N')
        seq_num_l = np.array([seq_dict[_] for _ in need_line_str]).astype(np.float32).transpose(1, 0)
        train_data_list.append(seq_num_l)
        train_label_list.append(np.array([int(data['label'])]))
    
    for idx, data in val_df.iterrows():
        need_line_str = str(data['sequence']).replace('\n', '')
        if len(need_line_str) <= seq_len:
            need_line_str = need_line_str.rjust(seq_len, 'N')
        seq_num_l = np.array([seq_dict[_] for _ in need_line_str]).astype(np.float32).transpose(1, 0)
        val_data_list.append(seq_num_l)
        val_label_list.append(np.array([int(data['label'])]))

    # Prepare DataLoader
    train_total_data = np.array(train_data_list)
    train_total_label = np.array(train_label_list)
    val_total_data = np.array(val_data_list)
    val_total_label = np.array(val_label_list)

    train_set = TensorDataset(torch.from_numpy(train_total_data).float(), torch.from_numpy(train_total_label).long())
    train_loader = DataLoader(dataset=train_set, batch_size=batch_size, shuffle=True, drop_last=True)
    val_set = TensorDataset(torch.from_numpy(val_total_data).float(), torch.from_numpy(val_total_label).long())
    val_loader = DataLoader(dataset=val_set, batch_size=batch_size, shuffle=True, drop_last=False)

    # Optimizer
    opt = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)

    # Logs for each epoch
    train_detail_dict = {'Epoch': [], 'Step': [], 'train_loss_l': [], 'train_acc_l': []}
    val_detail_dict = {'Epoch': [], 'Step': [], 'val_loss_l': [], 'val_acc_l': []}
    val_epoch_log_dict = {'Epoch': [], 'val_epoch_acc_l': []}

    # Training and validation loop
    best_val_acc = 0  # Initialize variable to track the best validation accuracy
    best_model = None  # Variable to store the best model

    for epoch in range(100):
        model.train()
        acc_num = 0
        loss_total = 0
        total_num = 0
        for train_data, train_label in train_loader:
            output = model(train_data)
            train_label = torch.squeeze(train_label)
            loss = F.cross_entropy(output, train_label)
            opt.zero_grad()
            loss.backward()
            opt.step()
            out_res = torch.argmax(output, dim=-1)
            acc_num = (out_res == train_label).sum().numpy()
            total_num = train_label.shape[0]
            loss_total = float(loss.detach().numpy())

        model.eval()
        acc_num = 0
        val_epoch_acc_num = 0
        for val_data, val_label in val_loader:
            output = model(val_data)
            out_res = torch.argmax(output, dim=-1)
            val_label = torch.squeeze(val_label)
            acc_num = (out_res == val_label).sum().numpy()
            val_epoch_acc_num += acc_num

        val_epoch_acc = val_epoch_acc_num / len(val_set)

        print(f'Epoch {epoch + 1}, Validation Accuracy: {val_epoch_acc * 100:.3f}%')

        # Save the model at the end of every epoch
        save_path = os.path.join(output_path, f'epoch_{epoch}.pkl')
        torch.save(model.state_dict(), save_path)

        # Track and save the best model based on validation accuracy
        if val_epoch_acc > best_val_acc:
            best_val_acc = val_epoch_acc
            best_model = model
            best_model_path = os.path.join(output_path, 'best_model.pkl')
            torch.save(best_model.state_dict(), best_model_path)  # Save the best model

        # Save epoch logs
        val_epoch_log_dict['Epoch'].append(epoch)
        val_epoch_log_dict['val_epoch_acc_l'].append(val_epoch_acc)

        # Save log files as needed
        val_epoch_df = pd.DataFrame(val_epoch_log_dict)
        val_epoch_df.to_csv(os.path.join(output_path, 'val_epoch_log.csv'), index=False)

    # After training, load the best model for final evaluation
    best_model = FCNA(motiflen=13)  # Reinitialize the model architecture
    best_model.load_state_dict(torch.load(best_model_path))

    # Evaluation: use the best model to evaluate performance on the validation set
    pre_seq = []
    labe_seq = []
    for val_data, val_label in val_loader:
        output = best_model(val_data)
        out_res = torch.argmax(output, dim=-1)
        pre_seq.append(out_res)
        labe_seq.append(torch.squeeze(val_label))

    pre_list = torch.cat(pre_seq).numpy()
    labe_seq = torch.cat(labe_seq).numpy()

    acc = (labe_seq == pre_list).sum() / labe_seq.shape[0]
    plot_confusion_matrix(labe_seq, pre_list, classes=[0, 1], title=f"PTK7, acc = {acc * 100:.3f}%")
    plt.savefig(f'{output_path}/2cls.pdf', bbox_inches='tight')
    print(f"Final Accuracy: {acc * 100:.3f}%")
    

    # Evaluate external data
    external_data_df = pd.read_csv(external_data, index_col=0)
    external_data_list, external_label_list = [], []

    for idx, data in external_data_df.iterrows():
        need_line_str = str(data['sequence']).replace('\n', '')
        if len(need_line_str) <= seq_len:
            need_line_str = need_line_str.rjust(seq_len, 'N')
        seq_num_l = np.array([seq_dict[_] for _ in need_line_str]).astype(np.float32).transpose(1, 0)
        external_data_list.append(seq_num_l)
        external_label_list.append(np.array([int(data['label'])]))

    external_total_data = np.array(external_data_list)
    external_total_label = np.array(external_label_list)
    external_set = TensorDataset(torch.from_numpy(external_total_data).float(), torch.from_numpy(external_total_label).long())
    external_loader = DataLoader(dataset=external_set, batch_size=batch_size, shuffle=False, drop_last=False)

    # Prediction and evaluation on external data
    pre_seq = []
    labe_seq = []
    for external_data, external_label in external_loader:
        output = best_model(external_data)
        out_res = torch.argmax(output, dim=-1)
        pre_seq.append(out_res)
        labe_seq.append(torch.squeeze(external_label))

    pre_list = torch.cat(pre_seq).numpy()
    labe_seq = torch.cat(labe_seq).numpy()

    acc = (labe_seq == pre_list).sum() / labe_seq.shape[0]

    res_df = pd.DataFrame({
        'True Label': labe_seq,
        'Predicted Label': pre_list
    })
    res_df.to_csv(f'{output_path}/external_validation_results.csv', index=False)   
    
    plot_confusion_matrix(labe_seq, pre_list, classes=[0, 1], title=f"PTK7 from Articles and Experiments\nacc = {acc * 100:.3f}%")
    plt.savefig(f"{output_path}/external_validation_confusion_matrix.pdf", bbox_inches="tight")
    print(f"External Data Accuracy: {acc * 100:.3f}%")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="FCNA Training and Evaluation")
    parser.add_argument('-train_data', type=str, required=True, help="Path to the label CSV file")
    parser.add_argument('-output_path', type=str, required=True, help="Directory path to save logs and models")
    parser.add_argument('-external_data', type=str, required=True, help="Path to the external data CSV file for evaluation")
    args = parser.parse_args()
    main(args)

