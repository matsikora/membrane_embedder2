from setuptools import setup
def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='membrane-embedder2',
      version='0.1',
      description='Create multiprotein systems for complexes',
      url='https://github.com/matsikora/membrane_embedder2',
      author='Mateusz Sikora',
      author_email='matsikora@gmail.com',
      license=' GNU GPLv3 ',
      packages=['membrane_embedder2'],
      install_requires=['numpy','mdanalysis','pytest'],
      zip_safe=False)
