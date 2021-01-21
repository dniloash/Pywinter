from numpy.distutils.core import setup, Extension

ext1 = Extension(name = 'KreatE_inter_m_f',
                 sources = ['pywinter/KreatE_inter_m_f.f90'])

setup(
  name = 'pywinter',
  packages = ['pywinter'],
  ext_modules = [ext1], # this must be the same as the name above
  version = '1.2.4',
  description = 'Create WRF-WPS intermediate files',
  long_description=("Python lib for creating WRF-WPS intermediate files.\n\n "
                      "Documentation:\n\n"
                      "https://pywinter.readthedocs.io/en/latest/\n"),
  author = 'Danilo A Suarez H',
  author_email = 'dniloash@gmail.com',
  url = 'https://github.com/dniloash/Pywinter',
  download_url = 'https://github.com/dniloash/Pywinter/tarball/1.2.4',
  keywords = ['Python3', 'WRF', 'WPS','Intermediate files'],
  classifiers = [],
  python_requires='>=3.6',
)
