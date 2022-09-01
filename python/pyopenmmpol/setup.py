import setuptools

setuptools.setup(
    name="pyopenmmpol",
    version="0.0.5",
    url="https://molimen1.dcci.unipi.it/filippo/open-mmpol/",
    author="Molecolab Group",
    author_email="benedetta.mennucci@unipi.it",
    description="Package that provide an interface to perform QM/MM calculations using atomistic polarizable embeddings",
    long_description=open('README.md').read(),
    packages=setuptools.find_packages(include=['pyopenmmpol']),
    install_requires=['numpy == 1.17.3'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
)

