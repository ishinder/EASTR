
import setuptools
setuptools.setup(
    name="EASTR",
    version="0.0.1",
    author="Ida Shinder",
    author_email="ishinde1@jhmi.edu",
    description="Tool for emending alignments of spuriously spliced transcript reads",
    url="https://github.com/ishinder/EASTR",
    install_requires = ["numpy","pandas","biopython","pysam","mappy"],
    packages=setuptools.find_packages(),
    entry_points={'console_scripts' : ['eastr = EASTR.run_eastr:main'], }
)