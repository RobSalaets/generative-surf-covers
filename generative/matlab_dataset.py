
import glob
from torch.utils.data import Dataset
from scipy.io import loadmat
import numpy
import torch


class MatlabDataset(Dataset):

    def __init__(self, keyword, transform=None):
        self.keyword = keyword
        self.transform = transform
        self.files = glob.glob("meshes/humans/*.mat")

    def __len__(self):
        return len(self.files)

    def __getitem__(self, idx):
        mat = loadmat(self.files[idx])[self.keyword]
        if self.transform:
            mat = self.transform(mat)
        return mat

class NumpyDataset(Dataset):

    def __init__(self):
        self.files = glob.glob("meshes/humans_numpy/*.npy")

    def __len__(self):
        return len(self.files)

    def __getitem__(self, idx):
        mat = numpy.load(self.files[idx])
        t = torch.from_numpy(mat)
        return t

