
import time
import random
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.data
import torchvision.datasets as dset
import torchvision.transforms as transforms
import torchvision.utils as vutils
import numpy as np
import matplotlib.pyplot as plt
from teeth_dataset import TeethDataset

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
        self.r1 = nn.LeakyReLU(negative_slope=0.2, inplace=True) # 4x4
        self.c2 = nn.ConvTranspose2d(hg, hg // 2, 5, 2, 2, 1, bias=False)
        self.bn2 = nn.BatchNorm2d(hg//2)
        # 8x8
        self.c3 = nn.ConvTranspose2d(hg // 2, hg // 4, 5, 2, 2, 1, bias=False)
        self.bn3 = nn.BatchNorm2d(hg // 4)
        # 16x16
        self.c4 = nn.ConvTranspose2d(hg // 4, hg // 8, 5, 2, 2, 1, bias=False)
        self.bn4 = nn.BatchNorm2d(hg // 8)
        # 32x32
        self.c5 = nn.ConvTranspose2d(hg // 8, nc, 5, 2, 2, 1, bias=False)
        self.tanh = nn.Tanh() # 64x64

    def forward(self, inp):
        oc1 = self.c1(inp)
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
            # input is (nc) x 64 x 64
            nn.Conv2d(nc, hd //8, 5, 2, 2, bias=False),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. (ndf) x 32 x 32
            nn.Conv2d(hd//8, hd//4, 5, 2, 2, bias=False),
            nn.BatchNorm2d(hd//4),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. (ndf*2) x 16 x 16
            nn.Conv2d(hd//4, hd//2, 5, 2, 2, bias=False),
            nn.BatchNorm2d(hd//2),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. (ndf*4) x 8 x 8
            nn.Conv2d(hd//2, hd, 5, 2, 2, bias=False),
            nn.BatchNorm2d(hd),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. (ndf*8) x 4 x 4
            nn.Conv2d(hd, 1, 4, 1, bias=False),
            nn.Sigmoid()
        )

    def forward(self, inp):
        return self.main(inp)


if __name__ == '__main__':
    manualSeed = random.randint(1, 10000)
    print("Random Seed: ", manualSeed)
    random.seed(manualSeed)
    torch.manual_seed(manualSeed)

    dataroot = "meshes/maps"
    workers = 2
    batch_size = 32
    image_size = 64
    orig_image_size = 100
    nc = 3 # Color channels
    nz = 100 # Latent vector size
    hg = 512 # number of feature maps
    hd = 512 # number of feature maps
    num_epochs = 100
    lr_g = 0.0001
    lr_d = 0.00007
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


    # netG = torch.load("netG1571059408.4237137.pt")
    netG = Generator(nz, image_size, nc, hg).to(device)
    netG.apply(weights_init)

    # netD = torch.load("netD1571059408.453633.pt")
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

    # Lists to keep track of progress
    img_list = []
    G_losses = []
    D_losses = []
    iters = 0

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
            if i % 50 == 0:
                print('[%d/%d][%d/%d]\tLoss_D: %.4f\tLoss_G: %.4f\tD(x): %.4f\tD(G(z)): %.4f / %.4f'
                      % (epoch, num_epochs, i, len(dataloader),
                         errD.item(), errG.item(), D_x, D_G_z1, D_G_z2))
            # if iters % 8 == 0:
            #     print("acc : %.4f" % last_acc)
            # Save Losses for plotting later
            G_losses.append(errG.item())
            D_losses.append(errD.item())

            # Check how the generator is doing by saving G's output on fixed_noise
            if (iters % 500 == 0) or ((epoch == num_epochs - 1) and (i == len(dataloader) - 1)):
                with torch.no_grad():
                    fake = netG(fixed_noise).detach().cpu()
                img_list.append(vutils.make_grid(fake, padding=2, normalize=True))
                # if i > 0:
                #     break

            iters += 1

    torch.save(netD, 'netD' + str(time.time() )+ '.pt')
    torch.save(netG, 'netG' + str(time.time() )+ '.pt')

    plt.figure(figsize=(8, 8))
    plt.axis("off")
    plt.imshow(
        np.transpose(vutils.make_grid(img_list[-1], padding=2, normalize=True).cpu(), (1, 2, 0)))
    plt.show()