from setuptools import setup, find_packages
import versioneer

__author__ = "Yan Hui"
__maintainer__ = "Yan Hui"
__copyright__ = "Copyright 2022, Yan Hui"
__email__ = "me@yanh.org"
__license__ = "GPL"

setup(
    name='laca',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url="https://github.com/yanhui09/laca",
    license=__license__,
    author=__author__,
    author_email=__email__,
    maintainer=__maintainer__,
    description="LACA, a reproducible and scalable workflow for Long Amplicon Consensus Analysis.",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
    ],
    entry_points={
        'console_scripts': [
            'laca = laca.laca:cli',
        ],
    },
     classifiers=["Topic :: Scientific/Engineering :: Bioinformatics"],
)