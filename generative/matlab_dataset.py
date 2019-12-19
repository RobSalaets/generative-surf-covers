
import glob
from torch.utils.data import Dataset
from scipy.io import loadmat


class MatlabDataset(Dataset):

    def __init__(self, size, channels, keyword, transform=None):
        self.size = size;
        self.channels = channels;
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
