import glob
from scipy.io import loadmat
import numpy
if __name__ == '__main__':
    files = glob.glob("meshes/humans/*.mat")
    for idx in range(len(files)):
        name = files[idx].split('\\')[1].split('.')[0] + '.npy'
        mat = loadmat(files[idx])['pushed_function'].astype(numpy.float32)
        mat = numpy.transpose(mat, [2, 0, 1])
        numpy.save('meshes/humans_numpy\\n' + name, mat)





