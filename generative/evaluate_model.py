import torch
import torch.nn as nn
import torch.utils.data
import torchvision.utils as vutils
import numpy as np
import matplotlib.pyplot as plt
from generative.train2 import Generator, Discriminator
from matplotlib.widgets import Slider
import scipy.io

import torchvision.transforms as transforms
from generative.matlab_dataset import MatlabDataset

if __name__ == '__main__':
    netG = torch.load("netG128_2.pt", map_location='cpu')
    netG.eval()
    netD = torch.load("netD128_2.pt", map_location='cpu')
    netD.eval()
    N = 64  #number of images
    dataset = MatlabDataset(128, 3, 'pushed_function',transform = transforms.Compose([transforms.ToTensor()]))
    dataloader = torch.utils.data.DataLoader(dataset, batch_size=32,shuffle=True, num_workers=1)
    device = torch.device("cuda:0")

    for i, data in enumerate(dataloader, 0):
    # # Format batch
        gpu_data = data.float().to(device)
    res2 = netD(gpu_data)

    fixed_noise = torch.randn(N, 128, 1, 1)
    fake = netG(fixed_noise)
    res = netD(fake)

    # for i, data in enumerate(dataloader, 0):
    #     # Format batch
    #     gpu_data = data.float().to(device)
    #     with torch.no_grad():
    #         # gpu_data = torch.nn.functional.interpolate(gpu_data, size=64, mode='bilinear', align_corners=False)
    #         plt.figure()
    #         plt.axis("off")
    #         plt.imshow(
    #             np.transpose(vutils.make_grid(gpu_data[0].detach().cpu(), padding=0, normalize=True).cpu(), (1, 2, 0)))
    #         plt.show()



    plt.figure(figsize=(4, 4))
    plt.axis("off")
    plt.imshow(
        np.transpose(vutils.make_grid(fake.detach(), padding=0, normalize=True).cpu(), (1, 2, 0)))
    plt.show()
    #
    #
    #
    # fig = plt.figure(figsize=(16, 8))
    # plt.axis("off")
    # im = plt.imshow(np.transpose(vutils.make_grid([fake[0].detach()], padding=2, normalize=True).cpu(), (1, 2, 0)))
    #
    # axcolor = 'lightgoldenrodyellow'
    # ax1 = plt.axes([0.25, 0.0, 0.65, 0.03], facecolor=axcolor)
    # ax2 = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=axcolor)
    # ax3 = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    # ax4 = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    # n1 = Slider(ax1, 'n1', 0.0, 1.0, valinit=.25, valstep=0.02)
    # n2 = Slider(ax2, 'n2', 0.0, 1.0, valinit=.25, valstep=0.02)
    # n3 = Slider(ax3, 'n3', 0.0, 1.0, valinit=.25, valstep=0.02)
    # n4 = Slider(ax4, 'n4', 0.0, 1.0, valinit=.25, valstep=0.02)


    # def update(val):
    #     sfn = (fixed_noise[:1] * n1.val + fixed_noise[1:2] * n2.val + fixed_noise[2:3] * n3.val + fixed_noise[
    #                                                                                               3:] * n4.val)
    #
    #     nfake = netG(sfn).detach()
    #     nearest = None
    #     im.set_array(np.transpose(vutils.make_grid(nfake[0], padding=2, normalize=True).cpu(), (1, 2, 0)))
    #     fig.canvas.draw_idle()


    # n1.on_changed(update)
    # n2.on_changed(update)
    # n3.on_changed(update)
    # n4.on_changed(update)
    # plt.show()