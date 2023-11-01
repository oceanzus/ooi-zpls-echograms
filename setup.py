from setuptools import setup, find_packages


# read the contents of the README file
def readme():
    with open('README.md', encoding='utf-8') as f:
        return f.read()


# read the contents of the VERSION file
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
      long_description_content_type='text/markdown',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Topic :: Data Parsing :: Scientific :: OOI :: Ocean Sonar :: Zooplankton',
      ],
      keywords=(
          'OOI Cabled Endurance Global Pioneer zooplankton echogram acoustic sonar ocean'
      ),
      url='https://github.com/oceanobservatories/ooi-zpls-echograms',
      author='Ocean Observatories Initiative',
      author_email='helpdesk@oceanobservatories.org',
      license='MIT',
      install_requires=[
          'bottleneck',
          'cmocean',
          'python-dateutil',
          'echopype>=0.8.1',
          'h5netcdf',
          'matplotlib',
          'numpy',
          'pandas',
          'pillow',
          'python-dateutil',
          'tqdm',
          'xarray'
      ],
      packages=find_packages(where="ooi_zpls_echograms"),
      package_dir={'': 'ooi_zpls_echograms'},
      include_package_data=True,
      zip_safe=False
      )
