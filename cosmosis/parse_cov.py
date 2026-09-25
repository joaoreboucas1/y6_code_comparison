import fitsio

filename = "EXTENSION_DVv8.0_fixmnu.fits"

with fitsio.FITS(filename, "r") as hdul:
    # Usually extension 1, but check with hdul if unsure
    print(hdul)

    covmat = hdul[1].read(columns=["COVMAT"])

print(covmat.shape)
output_file = "DESY6.cov"
with open(output_file, "w") as f:
    for i in range(covmat.shape[0]):
        for j in range(i, covmat.shape[1]):
            f.write(f"{i} {j} {covmat[i, j]}\n")
