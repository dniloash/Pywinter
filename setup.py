from numpy.distutils.core import setup, Extension

ext1 = Extension(name = 'KreatE_inter_m_f',
                 sources = ['pywinter/KreatE_inter_m_f.f90'])

setup(
  name = 'pywinter',
  packages = ['pywinter'],
  ext_modules = [ext1], # this must be the same as the name above
  version = '2.0.5',
  description = 'Read and Create WRF-WPS intermediate files',
  long_description=("Python lib for reading/creating WRF-WPS intermediate files.\n\n "
                      "Documentation:\n\n"
                      "https://pywinter.readthedocs.io/en/latest/\n"),
  author = 'Danilo A Suarez H',
  author_email = 'dniloash@gmail.com',
  url = 'https://github.com/dniloash/Pywinter',
  download_url = 'https://github.com/dniloash/Pywinter/tarball/2.0.5',
  keywords = ['Python', 'WRF', 'WPS','Intermediate files'],
  classifiers = [],
  python_requires='>=3.6',
)
