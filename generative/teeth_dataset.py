
import torch
import torch.nn
import glob
from torch.utils.data import Dataset
from scipy.io import loadmat
import numpy as np


class TeethDataset(Dataset):
    """Face Landmarks dataset."""

    def __init__(self, size, channels, root_dir, transform=None):
        self.size = size;
        self.channels = channels;
        self.root_dir = root_dir
        self.transform = transform
        self.files = glob.glob("meshes/maps/*p.mat")

    def __len__(self):
        return len(self.files)

    def __getitem__(self, idx):
        # if torch.is_tensor(idx):
        #     idx = idx.tolist()
        # if not isinstance(idx, list):
        #     idx = [idx];
        #
        # batch = np.empty((len(idx), self.channels, self.size, self.size), dtype=np.float32)
        # for i, ii in enumerate(idx):
        #     mat = loadmat(self.files[ii])['pushed_function']
        #     if self.transform:
        #         mat = self.transform(mat)
        #     batch[i, :, :, :] = mat
        mat = loadmat(self.files[idx])['tiled_rot']
        if self.transform:
            mat = self.transform(mat)
        return mat
