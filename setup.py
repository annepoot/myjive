from setuptools import setup, find_packages

setup(
    name="myjive",
    packages=find_packages(include=["myjive*", "myjivex*"]),
    version="0.1.8",
    license="MIT",
    description="Personal implementation of jive C++ library in Python",
    author="Anne Poot",
    author_email="a.poot-1@tudelft.nl",
    url="https://gitlab.tudelft.nl/apoot1/myjive",
    download_url="https://gitlab.tudelft.nl/apoot1/myjive/-/archive/v0.1.8/myjive-v0.1.8.tar.gz",
    keywords=[],
    python_requires="==3.10.*",
    install_requires=[
        "matplotlib==3.5.2",
        "numba==0.56.4",
        "numpy==1.23.5",
        "pandas==2.2.2",
        "scikit-sparse==0.4.12",
        "scipy==1.8.1",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
    ],
)
