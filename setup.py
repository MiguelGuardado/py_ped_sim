from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='ped_slim',
      version='0.1',
      description='ped_slim: a variational autoencoder for population genetic data',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/MiguelGuardado/ped_sim',
      author='Miguel Guardado',
      author_email='Miguel.Guardado@ucsf.edu',
      license='none',
      packages=find_packages(exclude=[]),
      scripts=["run_ped_sim.py"],
      zip_safe=False,
      )