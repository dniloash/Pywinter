from numpy.distutils.core import setup, Extension

ext1 = Extension(name = 'create_inter',
                 sources = ['pywinter/create_inter.f90'])

setup(
  name = 'pywinter',
  packages = ['pywinter'],
  ext_modules = [ext1], # this must be the same as the name above
  version = '1.1.0',
  description = 'Create WRF-WPS intermediate files',
  author = 'Danilo A Suarez H',
  author_email = 'dniloash@gmail.com',
  url = 'https://github.com/dniloash/Pywinter', 
  download_url = 'https://github.com/dniloash/Pywinter/tarball/1.1.0',
  keywords = ['Python3', 'WRF', 'WPS','Intermediate files'],
  classifiers = [],
  python_requires='>=3.6',
)
