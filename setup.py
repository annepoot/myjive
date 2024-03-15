from setuptools import setup, find_packages

setup(
  name = 'myjive',
  packages = find_packages(include=["myjive*", "myjivex*"]),
  version = '<X.X.X>',
  license='MIT',
  description = 'Personal implementation of jive C++ library in Python',
  author = 'Anne Poot',
  author_email = 'a.poot-1@tudelft.nl',
  url = 'https://gitlab.tudelft.nl/apoot1/myjive',
  download_url = 'https://gitlab.tudelft.nl/apoot1/myjive/-/archive/v<X.X.X>/myjive-v<X.X.X>.tar.gz',
  keywords = [],
  install_requires=[
    'matplotlib==3.5.2',
    'numba==0.56.4',
    'numpy==1.23.5',
    'pytest==7.1.2',
    'scikit-sparse==0.4.12',
    'scipy==1.8.1',
  ],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'Topic :: Scientific/Engineering :: Mathematics',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.10',
  ],
)
