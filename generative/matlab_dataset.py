
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
        self.files = glob.glob("meshes/humans/*.mat")

    def __len__(self):
        return len(self.files)

    def __getitem__(self, idx):
        mat = loadmat(self.files[idx])['tiled_rot']
        if self.transform:
            mat = self.transform(mat)
        return mat
