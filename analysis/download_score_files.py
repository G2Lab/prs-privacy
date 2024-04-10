from ftplib import FTP
import gzip
import os

# FTP server details
ftp_url = "ftp.ebi.ac.uk"
ftp_path = "/pub/databases/spot/pgs/scores/"
local_download_path = "catalog/"

# Connect to the FTP server
ftp = FTP(ftp_url)
ftp.login()

# List all folders on the server
# folders = ftp.nlst(ftp_path)
folders = [f"PGS00{i}" for i in range(4254, 4859)]
# print(folders)

# Iterate through folders
for folder in folders:
    folder_name = os.path.basename(folder)
    if "PGS" not in folder_name:
        continue
    print(f"Downloading {folder_name}")
    remote_file_path = f"{ftp_path}{folder}/ScoringFiles/Harmonized/{folder_name}_hmPOS_GRCh38.txt.gz"

    # Download the gzipped file
    local_file_path_gz = os.path.join(local_download_path, f"{folder_name}_hmPOS_GRCh38.txt.gz")
    with open(local_file_path_gz, "wb") as local_file:
        try:
            ftp.retrbinary(f"RETR {remote_file_path}", local_file.write)
        except Exception as e:
            print(f"{folder}: {e}")
            os.remove(local_file_path_gz)
            continue

    # Decompress the gzipped file
    with gzip.open(local_file_path_gz, "rb") as gz_file:
        file_content = gz_file.read()

    # Save the decompressed content to a new file
    local_file_path = os.path.join(local_download_path, f"{folder_name}_hmPOS_GRCh38.txt")
    with open(local_file_path, "wb") as local_file:
        local_file.write(file_content)

    # Remove the original gzipped file
    os.remove(local_file_path_gz)

# Close the FTP connection
ftp.quit()

# import aioftp
# import asyncio
# import os
# import gzip
#
#
# ftp_url = "ftp.ebi.ac.uk"
# ftp_path = "/pub/databases/spot/pgs/scores/"
# local_download_path = "catalog/"
#
#
# async def download_and_extract(sem, folder):
#     async with sem:
#         folder_name = os.path.basename(folder)
#         remote_file_path = f"{ftp_path}{folder}/ScoringFiles/Harmonized/{folder_name}_hmPOS_GRCh38.txt.gz"
#         print(f"Downloading {folder_name}")
#
#         # Connect to the FTP server
#         client = aioftp.Client()
#         await client.connect(ftp_url)
#
#         # Login to the FTP server (you might need credentials)
#         await client.login()
#
#         # Download the gzipped file
#         local_file_path_gz = os.path.join(local_download_path, f"{folder_name}_hmPOS_GRCh38.txt.gz")
#         await client.download(remote_file_path, local_file_path_gz)
#
#         # Disconnect from the FTP server
#         await client.quit()
#
#         # Decompress the gzipped file
#         with gzip.open(local_file_path_gz, "rb") as gz_file:
#             file_content = gz_file.read()
#
#         # Save the decompressed content to a new file
#         local_file_path = os.path.join(local_download_path, f"{folder_name}_hmPOS_GRCh38.txt")
#         with open(local_file_path, "wb") as local_file:
#             local_file.write(file_content)
#
#         # Remove the original gzipped file
#         os.remove(local_file_path_gz)
#
#
# async def main():
#     max_concurrent_tasks = 10
#     # List all folders on the server
#     # async with aioftp.Client.context("ftp.ebi.ac.uk") as client:
#     #     await client.login()
#     #     folders = await client.list("/pub/databases/spot/pgs/scores/")
#     folders = [f"PGS000{i}" for i in range(749, 1000)] + [f"PGS00{i}" for i in range(1000, 3979)]
#
#     # Create a semaphore to limit concurrent tasks
#     sem = asyncio.Semaphore(max_concurrent_tasks)
#
#     # Download and extract files concurrently
#     tasks = [download_and_extract(sem, folder) for folder in folders]
#     await asyncio.gather(*tasks)
#
# if __name__ == "__main__":
#     asyncio.run(main())

