from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()


def version():
    with open('VERSION.txt') as f:
        return f.read().strip()


setup(name='ooi_zpls_echograms',
      version=version(),
      description=(
          'Convert the raw data files from the OOI ocean sonar systems into '
          'processed NetCDF files and echogram plots (using the open source '
          'echopype package) to facilitate user access to, and exploration '
          'of, the bio-acoustic data collected across the multiple OOI arrays.'
      ),
      long_description=readme(),
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.8',
          'Topic :: Data Parsing :: Scientific :: OOI :: Ocean Sonar :: Zooplankton',
      ],
      keywords=(
          'OOI Cabled Endurance Global Pioneer zooplankton echogram acoustic sonar ocean'
      ),
      url='https://github.com/oceanobservatories/ooi-zpls-echograms',
      author='Ocean Observatories Initiative',
      author_email='helpdesk@oceanobservatories.org',
      license='MIT',
      packages=find_packages(),
      install_requires=[
          'cmocean',
          'python-dateutil',
          'echopype',
          'matplotlib',
          'numpy',
          'pandas',
          'pillow',
          'tqdm',
          'xarray'
      ],
      include_package_data=True,
      zip_safe=False)
