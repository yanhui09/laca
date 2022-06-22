from setuptools import setup, find_packages
import versioneer

__author__ = "Yan Hui"
__maintainer__ = "Yan Hui"
__copyright__ = "Copyright 2022, Yan Hui"
__email__ = "me@yanh.org"
__license__ = "GPL"

setup(
    name='kamp',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url="https://github.com/yanhui09/Kamp",
    license=__license__,
    author=__author__,
    author_email=__email__,
    maintainer=__maintainer__,
    description="Kamp, a k-mer based denoise pipeline for long read amplicon sequencing",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
    ],
    entry_points={
        'console_scripts': [
            'kamp = kamp.kamp:cli',
        ],
    },
     classifiers=["Topic :: Scientific/Engineering :: Bioinformatics"],
)