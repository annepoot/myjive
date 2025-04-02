from setuptools import setup, find_packages

setup(
    name="myjive",
    packages=find_packages(include=["myjive*", "myjivex*"]),
    version="0.1.9",
    license="MIT",
    description="Personal implementation of jive C++ library in Python",
    author="Anne Poot",
    author_email="a.poot-1@tudelft.nl",
    url="https://gitlab.tudelft.nl/apoot1/myjive",
    download_url="https://gitlab.tudelft.nl/apoot1/myjive/-/archive/v0.1.9/myjive-v0.1.9.tar.gz",
    keywords=[],
    python_requires="==3.10.*",
    install_requires=[
        "matplotlib",
        "numba",
        "numpy",
        "pandas",
        "scikit-sparse>=0.4.12",
        "scipy>=1.12",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
    ],
)
