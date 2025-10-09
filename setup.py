import os
from setuptools import setup, find_packages

# Function to read requirements from requirements.txt
def get_requires():
    reqs = []
    for line in open("CompactObject_TOV.egg-info/requires.txt", "r").readlines():
        reqs.append(line)
    return reqs

setup(
    name="CompactObject-TOV",
    version="2.1",
    packages=find_packages(),
    # install_requires=get_requires(),  # Automatically adds requirements
    # other setup parameters
)
