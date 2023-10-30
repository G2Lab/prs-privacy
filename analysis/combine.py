import fnmatch
import os


def combine_files(pgs_id, extension, dir):
    filepaths = [os.path.join(dir, filename) for filename in os.listdir(dir) if
                 fnmatch.fnmatch(filename, f"*{pgs_id}-*.{extension}") and "maf" not in filename and "eaf" not in filename]
    filepaths = sorted(filepaths, key=lambda x: int(x.split('-')[1].split(f".{extension}")[0]))
    print(filepaths)
    with open(os.path.join(dir, f"{pgs_id}.{extension}"), 'w') as outfile:
        for fname in filepaths:
            with open(fname) as infile:
                outfile.write(infile.read())
                # for line in infile:
                #     if line.startswith("individual"):
                #         continue
                #     outfile.write(line)


if __name__ == "__main__":
    combine_files("PGS002302", "json", "results/accuracyLikelihood")
