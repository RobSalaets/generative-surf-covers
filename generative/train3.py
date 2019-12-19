
import time
import random
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.data
import torchvision.datasets as dset
import torchvision.transforms as transforms
import torchvision.utils as vutils
import wandb
import numpy as np
import matplotlib.pyplot as plt
from teeth_dataset import TeethDataset
from layers import MinibatchStdDev, PeriodicConvTranspose2D, PeriodicConv2D

def weights_init(m):
    if type(m) == nn.Conv2d or type(m) == nn.ConvTranspose2d:
        nn.init.kaiming_normal_(m.weight)
    elif type(m) == nn.BatchNorm2d:
        nn.init.normal_(m.weight.data, 1.0, 0.02)
        nn.init.constant_(m.bias.data, 0)


class Generator(nn.Module):
    def __init__(self, nz, img_size, nc, hg):
        super(Generator, self).__init__()
        self.nz = nz
        self.img_size = img_size
        self.nc = nc
        self.hg = hg
        self.c1 = nn.ConvTranspose2d(nz, hg, 4, 1, 0, 0, bias=False)
        self.bn1 = nn.BatchNorm2d(hg)
        self.r1 = nn.LeakyReLU(negative_slope=0.2, inplace=True) # 4x4
        self.c2 = PeriodicConvTranspose2D(hg, hg // 2, 3, 2, 1, False)
        self.bn2 = nn.BatchNorm2d(hg//2)
        # 8x8
        self.c3 = PeriodicConvTranspose2D(hg // 2, hg // 4, 3, 2, 1, False)
        self.bn3 = nn.BatchNorm2d(hg // 4)
        self.c4 = PeriodicConvTranspose2D(hg // 4, hg // 8, 3, 2, 1, False)
        self.bn4 = nn.BatchNorm2d(hg // 8)

        # 16x16
        self.c5 = PeriodicConvTranspose2D(hg // 8, nc, 3, 2, 1, False)
        self.tanh = nn.Tanh() # 64x64

    def forward(self, inp):
        oc1 = self.c1(inp)
        oc1 = self.bn1(oc1)
        oc1 = self.r1(oc1)
        oc2 = self.c2(oc1)
        oc2 = self.bn2(oc2)
        oc2 = self.r1(oc2)
        oc3 = self.c3(oc2)
        oc3 = self.bn3(oc3)
        oc3 = self.r1(oc3)
        oc4 = self.c4(oc3)
        oc4 = self.bn4(oc4)
        oc4 = self.r1(oc4)
        oc5 = self.c5(oc4)

        return self.tanh(oc5)


class Discriminator(nn.Module):
    def __init__(self, img_size, nc, hd):
        super(Discriminator, self).__init__()
        self.img_size = img_size
        self.nc = nc
        self.hd = hd
        self.main = nn.Sequential(
            PeriodicConv2D(nc, hd//8, 3, 2, False),
            nn.LeakyReLU(0.2, inplace=True),
            PeriodicConv2D(hd // 8, hd//4, 3, 2, False),
            nn.BatchNorm2d(hd//4),
            nn.LeakyReLU(0.2, inplace=True),
            PeriodicConv2D(hd//4, hd//2, 3, 2, False),
            nn.BatchNorm2d(hd//2),
            nn.LeakyReLU(0.2, inplace=True),
            PeriodicConv2D(hd//2, hd, 3, 2, False),
            MinibatchStdDev(),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Conv2d(hd+1, 1, 4, 1, bias=False),
            nn.Sigmoid()
        )

    def forward(self, inp):
        return self.main(inp)


if __name__ == '__main__':
    manualSeed = random.randint(1, 10000)
    print("Random Seed: ", manualSeed)
    random.seed(manualSeed)
    torch.manual_seed(manualSeed)

    wandb.init(project="teeth-gan")

    dataroot = "meshes/maps"
    workers = 2
    batch_size = 20
    image_size = 64
    orig_image_size = 100
    nc = 3 # Color channels
    nz = 100 # Latent vector size
    hg = 128 # number of feature maps
    hd = 128 # number of feature maps
    num_epochs = 2000
    lr_sc = 2
    lr_g = 0.0005*lr_sc
    lr_d = 0.00001*lr_sc
    beta1 = 0.5 # Beta1 hyperparam for Adam optimizers

    dataset = TeethDataset(orig_image_size, nc, dataroot,
                           transform=transforms.Compose([
                               transforms.ToTensor()
                           ]))
    # Create the dataloader
    dataloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size,
                                             shuffle=True, num_workers=workers)
    device = torch.device("cuda:0")
    print(torch.cuda.get_device_name(torch.cuda.current_device()))


    # netG = torch.load("netG1576686127.5698137.pt")
    netG = Generator(nz, image_size, nc, hg).to(device)
    netG.apply(weights_init)
    #
    # netD = torch.load("netD1576686127.5539012.pt")
    netD = Discriminator(image_size, nc, hd).to(device)
    netD.apply(weights_init)

    criterion = nn.BCELoss()

    # Create batch of latent vectors that we will use to visualize
    #  the progression of the generator
    fixed_noise = torch.randn(64, nz, 1, 1, device=device)

    # Establish convention for real and fake labels during training
    real_label = .9
    fake_label = 0

    # Setup Adam optimizers for both G and D
    optimizerD = optim.Adam(netD.parameters(), lr=lr_d, betas=(beta1, 0.999))
    optimizerG = optim.Adam(netG.parameters(), lr=lr_g, betas=(beta1, 0.999))

    # Training Loop

    iters = 0

    wandb.watch(netG, log='all')
    wandb.watch(netD, log='all')

    print("Starting Training Loop...")
    # For each epoch
    for epoch in range(num_epochs):
        # For each batch in the dataloader
        last_acc = 0.0
        for i, data in enumerate(dataloader, 0):
            ############################
            # (1) Update D network: maximize log(D(x)) + log(1 - D(G(z)))
            ###########################
            ## Train with all-real batch
            netD.zero_grad()
            # Format batch
            gpu_data = data.float().to(device)
            with torch.no_grad():
                gpu_data = torch.nn.functional.interpolate(gpu_data, size=image_size, mode='bilinear', align_corners=False)
                gpu_data = gpu_data.mul(6.0)
            b_size = gpu_data.size(0)
            label = torch.full((b_size,), real_label, device=device)
            # Forward pass real batch through D
            output = netD(gpu_data).view(-1)
            # Calculate loss on all-real batch
            errD_real = criterion(output, label)
            # Calculate gradients for D in backward pass
            # if last_acc < .8:
            errD_real.backward()
            D_x = output.mean().item()

            ## Train with all-fake batch
            # Generate batch of latent vectors
            noise = torch.randn(b_size, nz, 1, 1, device=device)
            # Generate fake image batch with G
            fake = netG(noise)
            label.fill_(fake_label)
            # Classify all fake batch with D
            output = netD(fake.detach()).view(-1)
            # Calculate D's loss on the all-fake batch
            errD_fake = criterion(output, label)
            # Calculate the gradients for this batch
            # if last_acc < .8:
            errD_fake.backward()
            D_G_z1 = output.mean().item()
            # Add the gradients from the all-real and all-fake batches
            errD = errD_real + errD_fake
            # Update D
            # if last_acc < .8:
            optimizerD.step()
            # last_acc = (min(1.0, D_x / real_label) + 1.0 - D_G_z1) / 2.0

            ############################
            # (2) Update G network: maximize log(D(G(z)))
            ###########################
            netG.zero_grad()
            label.fill_(real_label)  # fake labels are real for generator cost
            # Since we just updated D, perform another forward pass of all-fake batch through D
            output = netD(fake).view(-1)
            # Calculate G's loss based on this output
            errG = criterion(output, label)
            # Calculate gradients for G
            errG.backward()
            D_G_z2 = output.mean().item()
            # Update G
            optimizerG.step()

            # Output training stats
            if i == 0:
                print('[%d/%d][%d/%d]\tLoss_D: %.4f\tLoss_G: %.4f\tD(x): %.4f\tD(G(z)): %.4f / %.4f'
                      % (epoch, num_epochs, i, len(dataloader),
                         errD.item(), errG.item(), D_x, D_G_z1, D_G_z2))
                wandb.log({'loss_D':errD.item(), 'loss_G':errG.item(), 'D(x)':D_x, 'D(G(z))_pre': D_G_z1, 'D(G(z))_post' : D_G_z2}, step=epoch)

            if iters % 100 == 0:
                plt.figure(figsize=(2,2))
                plt.axis("off")
                with torch.no_grad():
                    netG.eval()
                    fake = netG(fixed_noise).detach().cpu()
                    plt.imshow(
                        np.transpose(vutils.make_grid(fake[0:4], nrow=2,padding=2, normalize=True).cpu(), (1, 2, 0)))
                plt.show()
                # wandb.log({'examples': [wandb.Image(np.ones((64,64,3)), caption='ep%i'%(epoch))]}, step=epoch)
                netG.train()
            iters += 1


    netD._forward_hooks.clear()
    netG._forward_hooks.clear()
    torch.save(netD, 'netD' + str(time.time() )+ '.pt')
    torch.save(netG, 'netG' + str(time.time() )+ '.pt')

