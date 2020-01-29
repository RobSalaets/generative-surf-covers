import torchvision.transforms as transforms
from torch.utils.data import DataLoader
from generative.matlab_dataset import MatlabDataset, NumpyDataset
import torch.nn
from math import floor

class MatlabDataLoader:
    def __init__(self):
        # self.batch_table = {4:32, 8:32, 16:32, 32:16, 64:16, 128:16, 256:12, 512:3, 1024:1} # change this according to available gpu memory.
        self.batch_table = {4:64, 8:64, 16:64, 32:64, 64:64, 128:64, 256:64, 512:64, 1024:64} # change this according to available gpu memory.
        self.batchsize = int(self.batch_table[pow(2,2)])        # we start from 2^2=4
        self.imsize = int(pow(2,2))
        self.num_workers = 0
        # self.dataset = MatlabDataset('pushed_function',
        #                              transform=transforms.Compose([
        #                                  transforms.ToTensor()
        #                              ]))
        self.dataset = NumpyDataset()

    def renew(self, resl):
        print('[*] Renew dataloader configuration.')

        self.batchsize = int(self.batch_table[pow(2,resl)])
        self.imsize = int(pow(2,resl))


        self.dataloader = DataLoader(
            dataset=self.dataset,
            batch_size=self.batchsize,
            shuffle=True,
            num_workers=self.num_workers
        )

    def __iter__(self):
        return iter(self.dataloader)

    def __next__(self):
        return next(self.dataloader)

    def __len__(self):
        return len(self.dataloader.dataset)


    def get_batch(self, use_cuda):
        dataIter = iter(self.dataloader)
        batch = next(dataIter).float()
        if use_cuda:
            batch.cuda()
        return torch.nn.functional.interpolate(batch, size=self.imsize, mode='nearest')


        








