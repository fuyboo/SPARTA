
import logging
from enum import IntEnum, Enum
import click 
import numpy as np
from pathlib import Path

import torch
from torch import optim

from torch import nn
import torch.nn.functional as F
from tqdm import tqdm
import logging


def read_fasta(path):
    return get_reads_with_id_prefix(Path(path), ">", ">")


def read_fastq(path):
    return get_reads_with_id_prefix(Path(path), "@", "+")


def get_reads_with_id_prefix(path, prefix_on, prefix_off):
    reads = []
    read = ""
    switch = False
    with path.open() as f:
        for line in f.readlines():
            if line[0] == prefix_off:
                switch = False
                if len(read)>0 and not 'N' in read:
                    reads.append(read)
                read = ""
            if switch:
                read = read + line.strip()
            if line[0] == prefix_on:
                switch = True
                read = ""
        # write last read line
        if len(read)>0 and not 'N' in read:
            reads.append(read)
    return reads

class State(IntEnum):
    M = 0
    I = 1
    D = 2

class Transition(IntEnum):
    M2M = 0
    M2I = 1
    M2D = 2
    I2M = 3
    I2I = 4
    D2M = 5
    D2D = 6

def kld_loss(mu, logvar):
    KLD = - 0.5 * torch.sum(1 + logvar - mu.pow(2) -
                            logvar.exp()) / mu.shape[0]
    return KLD

class VAE(nn.Module):
    def __init__(self, encoder, decoder, embed_size=10, hidden_size=32):
        super(VAE, self).__init__()

        self.encoder = encoder
        self.decoder = decoder

        self.h2mu = nn.Linear(hidden_size, embed_size)
        self.h2logvar = nn.Linear(hidden_size, embed_size)

    def reparameterize(self, mu, logvar, deterministic=False):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        z = mu + (std * eps if not deterministic else 0)
        return z

    def forward(self, input, deterministic=False):
        h = self.encoder(input)
        mu = self.h2mu(h)
        logvar = self.h2logvar(h)

        z = self.reparameterize(mu, logvar, deterministic)
        recon_param = self.decoder(z)
        return recon_param, mu, logvar
    
class Bottleneck(nn.Module):
    def __init__(self, init_dim=32, window_size=7):
        super(Bottleneck, self).__init__()
        assert window_size % 2 == 1, f"window size should be odd, given {window_size}"

        self.conv1 = nn.Conv1d(
            in_channels=init_dim,
            out_channels=init_dim*2,
            kernel_size=1)

        self.conv2 = nn.Conv1d(
            in_channels=init_dim*2,
            out_channels=init_dim*2,
            kernel_size=window_size,
            padding=window_size//2
        )

        self.conv3 = nn.Conv1d(
            in_channels=init_dim*2,
            out_channels=init_dim,
            kernel_size=1)

        self.bn1 = nn.BatchNorm1d(init_dim)
        self.bn2 = nn.BatchNorm1d(init_dim*2)
        self.bn3 = nn.BatchNorm1d(init_dim*2)

    def forward(self, input):
        x = self.conv1(F.leaky_relu(self.bn1(input)))
        x = self.conv2(F.leaky_relu(self.bn2(x)))
        x = self.conv3(F.leaky_relu(self.bn3(x)))
        return F.leaky_relu(x+input)
    
class DecoderPHMM(nn.Module):
    # tile hidden and input to make x
    def __init__(self,  motif_len, embed_size,  hidden_size=32):
        super(DecoderPHMM, self).__init__()

        class View(nn.Module):
            def __init__(self, shape):
                super(View, self).__init__()
                self.shape = shape

            def forward(self, x):
                return x.view(*self.shape)

        self.fc1 = nn.Sequential(
            nn.Linear(embed_size, hidden_size),
            nn.BatchNorm1d(hidden_size),
            nn.LeakyReLU(negative_slope=0.01, inplace=True))

        self.fc2 = nn.Sequential(
            nn.Linear(hidden_size, hidden_size*2),
            nn.BatchNorm1d(hidden_size*2),
            nn.LeakyReLU(negative_slope=0.01),
            nn.Linear(hidden_size*2, hidden_size),
            nn.BatchNorm1d(hidden_size),
            nn.LeakyReLU(negative_slope=0.01)
        )

        self.tr_from_M = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            nn.Linear(hidden_size, (motif_len+1)*3),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            View((-1, motif_len+1, 3)),
            nn.LogSoftmax(dim=2)
        )
        self.tr_from_I = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            nn.Linear(hidden_size, (motif_len+1)*2),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            View((-1, motif_len+1, 2)),
            nn.LogSoftmax(dim=2)
        )
        self.tr_from_D = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            nn.Linear(hidden_size, (motif_len+1)*2),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            View((-1, motif_len+1, 2)),
            nn.LogSoftmax(dim=2)
        )

        self.emission = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            nn.Linear(hidden_size, motif_len*4),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            View((-1, motif_len, 4)),
            nn.LogSoftmax(dim=2)
        )

    def forward(self, input):
        x = self.fc1(input)

        transition_from_match = self.tr_from_M(x)
        transition_from_insertion = self.tr_from_I(x)
        transition_from_deletion = self.tr_from_D(x)

        emission_proba = self.emission(x)
        return (torch.cat((
            transition_from_match,
            transition_from_insertion,
            transition_from_deletion), dim=2), emission_proba)


class EncoderCNN (nn.Module):
    # 0~3 is already used by embedding ATGC
    def __init__(self, embedding_dim=32, window_size=7, num_layers=6):
        super(EncoderCNN, self).__init__()
        self.embedding_dim = embedding_dim
        self.window_size = window_size

        self.embed = nn.Embedding(
            num_embeddings=4,  # [A,T,G,C,PAD,SOS,EOS]
            embedding_dim=embedding_dim)

        modules = [Bottleneck(embedding_dim, window_size)
                   for _ in range(num_layers)]
        self.resnet = nn.Sequential(*modules)

    def forward(self, seqences):
        # change X from (N, L) to (N, L, C)
        x = F.leaky_relu(self.embed(seqences))

        # change X to (N, C, L)
        x = x.transpose(1, 2)
        value, indices = self.resnet(x).max(dim=2)
        return value
    
class CNN_PHMM_VAE(VAE):
    def __init__(self, motif_len=12, embed_size=10, hidden_size=32, kernel_size=7):
        encoder = EncoderCNN(hidden_size, kernel_size)
        decoder = DecoderPHMM(motif_len, embed_size)

        super(CNN_PHMM_VAE, self).__init__(
            encoder, decoder, embed_size, hidden_size)
        self.loss_fn = profile_hmm_loss_fn

def profile_hmm_loss(recon_param, input, force_matching=False, match_cost=5):
    batch_size, random_len = input.shape

    a, e_m = recon_param
    motif_len = e_m.shape[1]

    F = torch.ones((batch_size, 3, motif_len + 1, random_len + 1),
                   device=input.device) * (-100)
    # init
    F[:, 0, 0, 0] = 0

    for i in range(random_len + 1):
        for j in range(motif_len + 1):
            # State M
            if j*i != 0:
                F[:, State.M, j, i] = e_m[:, j - 1].gather(1, input[:, i - 1:i])[:, 0] + \
                    torch.logsumexp(torch.stack((
                        a[:, j - 1, Transition.M2M] +
                        F[:, State.M, j - 1, i - 1],
                        a[:, j - 1, Transition.I2M] +
                        F[:, State.I, j - 1, i - 1],
                        a[:, j - 1, Transition.D2M] +
                        F[:, State.D, j - 1, i - 1])), dim=0)

            # State I
            if i != 0:
                F[:, State.I, j, i] = - 1.3863 + \
                    torch.logsumexp(torch.stack((
                        a[:, j, Transition.M2I] +
                        F[:, State.M, j, i-1],
                        # Removed D-to-I transition
                        # a[:, j, Transition.D2I] +
                        # F[:, State.D, j, i-1],
                        a[:, j, Transition.I2I] +
                        F[:, State.I, j, i-1]
                    )), dim=0)

            # State D
            if j != 0:
                F[:, State.D, j, i] = \
                    torch.logsumexp(torch.stack((
                        a[:, j - 1, Transition.M2D] +
                        F[:, State.M, j - 1, i],
                        # REMOVED I-to-D transition
                        # a[:, j - 1, Transition.I2D] +
                        # F[:, State.I, j - 1, i],
                        a[:, j - 1, Transition.D2D] +
                        F[:, State.D, j - 1, i]
                    )), dim=0)

    # final I->M transition
    F[:, State.M, motif_len, random_len] += a[:,
                                              motif_len, Transition.M2M]
    F[:, State.I, motif_len, random_len] += a[:,
                                              motif_len, Transition.I2M]
    F[:, State.D, motif_len, random_len] += a[:,
                                              motif_len, Transition.D2M]

    if force_matching:
        force_loss = np.log((match_cost+1)*match_cost/2) + \
            torch.sum((match_cost-1) * a[:, :, Transition.M2M], dim=1).mean()
        return - force_loss - torch.logsumexp(F[:, :, motif_len, random_len], dim=1).mean()
    return - torch.logsumexp(F[:, :, motif_len, random_len], dim=1).mean()


def profile_hmm_loss_fn(input, recon_param, mu, logvar, debug=False, test=False, beta=1, force_matching=False, match_cost=5):
    phmmloss = profile_hmm_loss(
        recon_param, input, force_matching=force_matching, match_cost=match_cost)
    kld = kld_loss(mu, logvar)

    if test:
        return phmmloss.item(), kld.item()
    return phmmloss + beta * kld
class nt_index(IntEnum):
    A = 0
    T = 1
    G = 2
    C = 3
    PAD = 4
    SOS = 5
    EOS = 6
    U = 1

def one_hot_index(seq):
    return [int(nt_index[char]) for char in seq]


def one_hot_encode(nucleotide, padding=0):
    arr = np.vstack((np.eye(4), np.ones(4)[None, :]*0.25))
    return arr[one_hot_index("N"*padding + nucleotide + "N"*padding)].T

def train(epochs, model, train_loader, test_loader, optimizer,logger, loss_fn=None, device="cuda", n_print=100, model_str="model.mdl", save_dir=Path("./"), threshold=20, beta=1, beta_schedule=False, force_matching=False, force_epochs=20, logs=True, position=0):
    csv_filename = model_str.replace(".mdl", ".csv")
    if loss_fn == profile_hmm_loss_fn and force_matching:
        logger.info(f"force till {force_epochs}")
    patient = 0
    losses = []
    test_losses = []

    if loss_fn is None:
        loss_fn = model.loss_fn

    with tqdm(total=epochs, position=position+1) as pbar:
        description = ""
        for epoch in range(1, epochs + 1):
            if beta_schedule and epoch < threshold:
                beta = epoch / threshold
            model.train()
            train_loss = 0
            for data in train_loader:
                data = data.to(device)
                optimizer.zero_grad()
                if loss_fn in {profile_hmm_loss_fn, profile_hmm_loss_fn_fast} and epoch <= force_epochs:
                    loss = loss_fn(data, *model(data), beta=beta,
                                   force_matching=force_matching, match_cost=1+4*(1-epoch/force_epochs))
                else:
                    loss = loss_fn(data, *model(data), beta=beta)
                loss.backward()
                train_loss += loss.item() * data.shape[0]
                optimizer.step()
            train_loss /= len(train_loader.dataset)
            if np.isnan(train_loss):
                logger.info("!-- train -> nan")
                if len(losses) == 0:
                    return [(float("inf"), float("inf"), float("inf"), float("inf"))]
                return losses
            model.eval()
            test_ce = 0
            test_kld = 0
            with torch.no_grad():
                for data in test_loader:
                    data = data.to(device)
                    ce, kld = loss_fn(data, *model(data), beta=beta, test=True)
                    test_ce += ce * data.shape[0]
                    test_kld += kld * data.shape[0]
            test_kld /= len(test_loader.dataset)
            test_ce /= len(test_loader.dataset)
            test_loss = test_kld + test_ce
            if np.isnan(test_loss):
                logger.info("!-- test -> nan")
                if len(losses) == 0:
                    return [(float("inf"), float("inf"), float("inf"), float("inf"))]
                return losses
            loss = (train_loss, test_loss, test_ce, test_kld)
            losses.append(loss)
            test_losses.append(test_loss)
            if len(test_losses) - 1 == np.argmin(test_losses):
                torch.save(model.state_dict(), save_dir / model_str)
                patient = 0
            else:
                patient += 1
                if patient > threshold:
                    # logger.info(f"{epoch}: no progress in test loss for {patient} iteration. break.")
                    return losses

            patience_str = f"[{patient}]" if patient > 0 else (
                "[" + "⠸⠴⠦⠇⠋⠙"[epoch % 6] + "]")
            len_model_str = len(model_str)
            if len_model_str > 10:
                model_str_print = f"..........{model_str}.........."[
                    (epoch+9) % (len_model_str+10):(epoch+9) % (len_model_str+10) + 10]
            else:
                model_str_print = model_str
            description = f'{patience_str:>4}{epoch:4d} itr {train_loss:6.2f} <-> {test_loss:6.2f} ({test_ce:6.2f}+{test_kld:6.2f}) of {model_str_print}'

            if epoch == 1:
                with open(save_dir / csv_filename, "w") as f:
                    f.write("epoch,train_loss,test_loss,test_recon,test_kld\n")
            with open(save_dir / csv_filename, "a") as f:
                f.write(f"{epoch}," + ",".join(map(str, loss)) + "\n")
            if logs:
                logger.debug(description)
                pbar.set_description(description)
                pbar.update(1)
    return losses

def profile_hmm_loss_fn_fast(input, recon_param, mu, logvar, debug=False, test=False, beta=1, force_matching=False, match_cost=5):
    phmmloss = torch_multi_polytope_dp_log(*recon_param, input, force_matching, match_cost)
    kld = kld_loss(mu, logvar)

    if debug:
        logger.info(f"phmm={phmmloss:.2f}, kld={kld:.2f}")
    if test:
        return phmmloss.item(), kld.item()
    return phmmloss + beta * kld

def torch_multi_polytope_dp_log(transition_proba, emission_proba, output, force_matching=False, match_cost=5):
    """
    torch_multi_polytope_dp_log(
        transition_proba,
        emission_proba,
        output
    )

    Given logarithmic parameters, the function calculates 
    the probability of the output sequence of a certain 
    Profile Hidden Markov Model (PHMM) by forward algorithm.

    For the efficiency, this function is utilizing polytope 
    model which enables parallel dynamic programming (DP).

    Parameters
    ----------
    transition_proba : torch.Tensor
        the tensor which define the probability to transit
        state to state. the tensor shape has to be 
        (`batch`, `from`=3, `to`=3, `model_length`+1) and
        the tensor has to be logarithmic number

    emission_proba : torch.Tensor
        the tensor which emit characters. The tensor shape 
        has to be: (`batch`, `model_length`, `augc`=4)

    output : torch.Tensor
        the tensor of the output vector. the tensor shape
        has to be (`batch`, `string_length`)

    Returns
    ------
    probabilities : torch.Tensor
        log-probabilities of the given output tensor with
        shape (`batch`,)

    """
    model_length = emission_proba.shape[1]
    batch_size, string_length = output.shape

    F = torch.ones(
        size=(batch_size, 3, model_length +
              string_length + 1, string_length + 1),
        device=output.device) * - 200
    F[:, State.M, 0, 0] = 0
    log4 = torch.Tensor([4]).log().to(output.device)
    arange = torch.arange(
        start=0, end=model_length + string_length + 1, device=output.device)
    for model_index_pre in range(1, model_length + string_length+1):
        if max(1, model_index_pre-model_length) < min(string_length+1, model_index_pre):
            m_slice = arange[max(1, model_index_pre-model_length)                             : min(string_length+1, model_index_pre)]
            F[:, State.M, model_index_pre, m_slice] = \
                torch.gather(
                    emission_proba[:, model_index_pre - m_slice - 1], 2,
                    output[:, m_slice - 1, None]).reshape(batch_size, len(m_slice)) \
                + torch.logsumexp(
                    transition_proba[:, :, State.M,
                                     model_index_pre - m_slice - 1]
                    + F[:, :, model_index_pre - 2, m_slice - 1], axis=1)

        if max(1, model_index_pre-model_length) < min(string_length+1, model_index_pre+1):
            i_slice = arange[max(1, model_index_pre-model_length)                             : min(string_length+1, model_index_pre+1)]
            F[:, State.I, model_index_pre, i_slice] = \
                torch.logsumexp(
                    transition_proba[:, :, State.I, model_index_pre - i_slice]
                    + F[:, :, model_index_pre - 1, i_slice - 1], axis=1) \
                - log4

        if max(0, model_index_pre-model_length) < min(string_length+1, model_index_pre):
            d_slice = arange[max(1, model_index_pre-model_length)                             : min(string_length+1, model_index_pre)]
            F[:, State.D, model_index_pre, d_slice] = \
                torch.logsumexp(
                    transition_proba[:, :, State.D,
                                     model_index_pre - d_slice - 1]
                    + F[:, :, model_index_pre - 1, d_slice], axis=1)
    if force_matching:
        return -torch.logsumexp(F[:, :, -1, -1] + transition_proba[:, :, State.M, -1], axis=1).mean()\
            - np.log((match_cost+1)*match_cost/2) \
            - torch.sum((match_cost-1) * transition_proba[:, State.M, State.M, :], dim=1).mean()
    return -torch.logsumexp(F[:, :, -1, -1] + transition_proba[:, :, State.M, -1], axis=1).mean()

class ProfileHMMSampler():
    def __init__(self, transition_proba, emission_proba, proba_is_log=False):
        self.e = emission_proba
        self.a = transition_proba
        if proba_is_log:
            self.e = np.exp(self.e)
            self.a = np.exp(self.a)
        self.e = self.e / np.sum(self.e, axis=1)[:, None]

    def sample(self, sequence_only=False, debug=False):
        idx, state = (0, State.M)
        states = [(idx, state)]
        seq = ""
        while True:
            if state == State.M:
                p = self.a[idx][np.array([
                    Transition.M2M.value,
                    Transition.M2I.value,
                    Transition.M2D.value])]
            elif state == State.I:
                p = np.stack([
                    self.a[idx][Transition.I2M.value],
                    self.a[idx][Transition.I2I.value],
                    0])
            elif state == State.D:
                p = np.stack([
                    self.a[idx][Transition.D2M.value],
                    0,
                    self.a[idx][Transition.D2D.value]])
            else:
                logger.info("something wrong")

            state = np.random.choice([State.M, State.I, State.D], p=p/sum(p))
            if state != State.I:
                idx += 1
            states.append((idx, state))
            if idx == self.a.shape[0]:
                break

            if state == State.M:
                # logger.info("{:.2f}, {:.2f}, {:.2f}, {:.2f}".format(*self.e[idx-1]))

                seq += np.random.choice(list("ATGC"), p=self.e[idx-1])
                if debug:
                    logger.info(idx, state, self.e[idx-1], seq[-1])
            elif state == State.I:
                seq += np.random.choice(list("atgc"))
            else:
                seq += "_"
        if not sequence_only:
            return states, seq
        else:
            return seq

    def most_probable(self, sequence_only=False):
        idx, state = (0, State.M)
        states = [(idx, state)]
        seq = ""
        while True:
            if state == State.M:
                p = self.a[idx][np.array([
                    Transition.M2M.value,
                    Transition.M2I.value,
                    Transition.M2D.value])]
            elif state == State.I:
                p = [
                    self.a[idx][Transition.I2M.value],
                    0,
                    0]
            elif state == State.D:
                p = [
                    self.a[idx][Transition.D2M.value],
                    0,
                    self.a[idx][Transition.D2D.value]]
            else:
                logger.info("something wrong")
            p[np.argmax(p)] += 1000000
            state = np.random.choice([State.M, State.I, State.D], p=p/sum(p))
            if state != State.I:
                idx += 1
            states.append((idx, state))

            if idx == self.a.shape[0]:
                break

            if state == State.M:
                # logger.info("{:.2f}, {:.2f}, {:.2f}, {:.2f}".format(*self.e[idx-1]))
                p = np.copy(self.e[idx-1])
                p[np.argmax(p)] += 100000
                seq += np.random.choice(list("ATGC"), p=p/sum(p))
            elif state == State.I:
                seq += "N"
            else:
                seq += "_"
        if not sequence_only:
            return states, seq
        else:
            return seq

    def calc_seq_proba(self, seq: str):
        one_hot_seq = torch.tensor(one_hot_index(seq))
        model_len = self.e.shape[0]
        random_len = len(seq)

        e = np.log(self.e)
        a = np.log(self.a)

        F = torch.ones((3, model_len + 2, random_len + 1)) * (-100)

        # init
        F[0, 0, 0] = 0

        for i in range(random_len + 1):
            for j in range(model_len + 1):
                # State M
                if j*i != 0:
                    F[State.M, j, i] = e[j - 1][one_hot_seq[i - 1]] + \
                        torch.logsumexp(torch.stack((
                            a[j - 1, Transition.M2M] +
                            F[State.M, j - 1, i - 1],
                            a[j - 1, Transition.I2M] +
                            F[State.I, j - 1, i - 1],
                            a[j - 1, Transition.D2M] + F[State.D, j - 1, i - 1])), dim=0)

                # State I
                if i != 0:
                    F[State.I, j, i] = - 1.3863 + \
                        torch.logsumexp(torch.stack((
                            a[j, Transition.M2I] + F[State.M, j, i-1],
                            a[j, Transition.I2I] + F[State.I, j, i-1]
                        )), dim=0)

                # State D
                if j != 0:
                    F[State.D, j, i] = \
                        torch.logsumexp(torch.stack((
                            a[j - 1, Transition.M2D] + F[State.M, j - 1, i],
                            a[j - 1, Transition.D2D] + F[State.D, j - 1, i]
                        )), dim=0)

        F[State.M, model_len+1, random_len] = \
            torch.logsumexp(torch.stack((
                a[model_len, Transition.M2M] +
                F[State.M, model_len, random_len],
                a[model_len, Transition.I2M] +
                F[State.I, model_len, random_len],
                a[model_len, Transition.D2M] +
                F[State.D, model_len, random_len]
            )), dim=0)

        return F[State.M, model_len+1, random_len]