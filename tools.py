
from alphapept.pyrawfilereader import RawFileReader

fpath = '/home/mleeming/MSPF/MSPF_SHARES/MEDIAFLUX/EXPLORIS_3/210419_AlexA_DIA1_4.raw'
rawfile = RawFileReader(fpath)

spec_indices = range(rawfile.FirstSpectrumNumber, rawfile.LastSpectrumNumber + 1)


for ix, i in enumerate(spec_indices):

    rt_1 = rawfile.RTFromScanNum(i)
    order_1 = rawfile.GetMSOrderForScanNum(i)

    rt_2 = rawfile.RTFromScanNum(i+1)
    order_2 = rawfile.GetMSOrderForScanNum(i+1)

#    if 1 in [order_1, order_2]:
#        print( order_2, order_1, (rt_2 - rt_1) * 60)

    if order_1 == 2:
        mzs, ints = rawfile.GetProfileMassListFromScanNum(i)

        for j in range(1, len(mzs)):
            print(mzs[j] - mzs[j - 1])

