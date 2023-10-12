import os
import sys

def create_folder(folder_name):
    os.makedirs(folder_name, exist_ok=True)

if __name__ == "__main__":

    # Create a folder
    folder_name = "output_files"
    create_folder(folder_name)
