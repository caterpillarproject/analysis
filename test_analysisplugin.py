import haloutils
from analysisplugin import NvmaxPlugin,SHMFPlugin

if __name__=="__main__":
    Nvmax = NvmaxPlugin()
    SHMF = SHMFPlugin()
    for hpath in ["/bigbang/data/AnnaGroup/caterpillar/halos/H95289/H95289_BB_Z127_P7_LN7_LX11_O4_NV4_B",
                  "/bigbang/data/AnnaGroup/caterpillar/halos/H95289/H95289_BB_Z127_P7_LN7_LX12_O4_NV4_B"]:
        Nvmax(hpath)
        SHMF(hpath)
