from setuptools import setup, find_packages

setup(
    name='GSadjust',
    version='2.0',
    description='Gravity-survey network adjustment',
    author='',
    author_email='jkennedy@usgs.gov',
    url='https://github.com/jkennedy-usgs/sgp-gsadjust',
    packages=find_packages(exclude=['test_data', 'docs', 'tests*']),
)
