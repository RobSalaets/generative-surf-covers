import torch
import torch.nn as nn

class MinibatchStdDev(nn.Module):

    def __init__(self):
        super().__init__()

    def forward(self, x):
        """
        :param x: shape [N, C, H, W]
        """
        (N,C,H,W) = x.shape
        std = torch.mean(torch.std(x, 0), (0, 1, 2))
        map = torch.full((N, 1, H, W), std.item())
        if x.is_cuda:
            map = map.cuda()
        x = torch.cat((x,map), 1)
        return x

class PeriodicConvTranspose2D(nn.Module):
    def __init__(self, fin, fout, kernel, stride, output_padding, use_bias):
        super().__init__()
        self.fin = fin
        self.fout = fout
        self.kernel = kernel
        self.stride = stride
        self.output_padding = output_padding
        self.use_bias = use_bias
        assert(stride==2)
        assert(kernel == 5 or kernel == 3)
        self.left_pad = 0 if kernel == 3 else 1
        self.top_pad =  0 if kernel == 3 else 1
        self.right_pad = 1 if kernel == 3 else 1
        self.bottom_pad = 1 if kernel == 3 else 1
        self.half_kernel = (kernel-1)//2
        self.ct2d = nn.ConvTranspose2d(fin, fout, kernel, stride, 0, 0, bias=use_bias)
        # self.ct2d = nn.ConvTranspose2d(fin, fout, kernel, stride, kernel-1+(kernel-1)//2, output_padding, bias=use_bias)

    @property
    def weight(self):
        return self.ct2d.weight

    @weight.setter
    def weight(self, w):
        self.ct2d.weight = w

    def forward(self, x):
        """
        :param x: shape [N, C, H, W]
        """
        if self.top_pad == 1:
            slice_bottom = torch.cat((x[:,:,-self.top_pad:,-self.left_pad:], x[:,:,-self.top_pad:,:], x[:,:,-self.top_pad:, 0:self.right_pad]), 3)
            slice_top = torch.cat((x[:,:,0:self.bottom_pad,-self.left_pad:],x[:,:,0:self.bottom_pad,:],x[:,:,0:self.bottom_pad,0:self.right_pad]), 3)
            mid = torch.cat((x[:,:,:,-self.left_pad:],x, x[:,:,:,0:self.right_pad]), 3)
            padded = torch.cat((slice_bottom, mid, slice_top), 2)
        else:
            slice_top = torch.cat((x[:, :, 0:self.bottom_pad, :],
                                   x[:, :, 0:self.bottom_pad, 0:self.right_pad]), 3)
            mid = torch.cat((x, x[:, :, :, 0:self.right_pad]), 3)
            padded = torch.cat((mid, slice_top), 2)
        out = self.ct2d(padded)
        return out[:,:, self.half_kernel:-(self.half_kernel+1), self.half_kernel:-(self.half_kernel+1)]


class PeriodicConv2D(nn.Module):
    def __init__(self, fin, fout, kernel, stride, use_bias):
        super().__init__()
        self.fin = fin
        self.fout = fout
        self.kernel = kernel
        self.stride = stride
        self.use_bias = use_bias

        # Assumes stride 2
        assert(stride <=2)
        if stride == 2:
            self.left_pad = (kernel-1) // 2
            self.top_pad = (kernel-1) // 2
            self.right_pad = (kernel-1) // 2 -1
            self.bottom_pad = (kernel-1) // 2 -1
        else:
            self.left_pad = (kernel - 1) // 2
            self.top_pad = (kernel - 1) // 2
            self.right_pad = (kernel - 1) // 2
            self.bottom_pad = (kernel - 1) // 2
        self.c2d = nn.Conv2d(fin, fout, kernel, stride, 0, bias=use_bias)

    @property
    def weight(self):
        return self.c2d.weight

    @weight.setter
    def weight(self, w):
        self.c2d.weight = w

    def forward(self, x):
        """
        :param x: shape [N, C, H, W]
        """

        slice_bottom = torch.cat((x[:,:,-self.top_pad:,-self.left_pad:], x[:,:,-self.top_pad:,:], x[:,:,-self.top_pad:, 0:self.right_pad]), 3)
        slice_top = torch.cat((x[:,:,0:self.bottom_pad,-self.left_pad:],x[:,:,0:self.bottom_pad,:],x[:,:,0:self.bottom_pad,0:self.right_pad]), 3)
        mid = torch.cat((x[:,:,:,-self.left_pad:], x, x[:,:,:,0:self.right_pad]), 3)
        padded = torch.cat((slice_bottom, mid, slice_top), 2)
        out = self.c2d(padded)
        return out