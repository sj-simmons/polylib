from setuptools import setup, find_packages
#from mypyc.build import mypycify

def readme():
  with open('README.rst', 'r') as fh:
    return fh.read()

setup(
  entry_points={
      'console_scripts': [
          'bernoulli = polylib.bernoulli:main',
      ],
      'gui_scripts': [],
  },
  name='polylib',
  url='https://github.com/sj-simmons/polylib',
  download_url='https://github.com/sj-simmons/polylib/archive/v0.3.1.tar.gz',
  author='Scott Simmons',
  author_email='ssimmons@drury.edu',
  packages=find_packages(),
  #ext_modules=mypycify([
  #    'polylib/polynomials.py',
  #    #'polylib/bernoulli.py',
  #    ]),
  python_requires='>=3.8',
  install_requires=['numlib'],
  version="0.3.1",
  license='Apache 2.0',
  description='a library for working with polynomials',
  long_description=readme(),
  include_package_data=True,
  zip_safe=False,
  project_urls={'Upstream Repository': 'https://gihub.com/sj-simmons/polylib'},
  classifiers=[
      'Development Status :: 4 - Beta',
      'Intended Audience :: Science/Research',
      'Intended Audience :: Education',
      'Programming Language :: Python :: 3.8',
      'Programming Language :: Python :: 3.9',
      'Programming Language :: Python :: 3.10',
      'Programming Language :: Python :: 3.11',
      'Topic :: Scientific/Engineering',
      'Topic :: Scientific/Engineering :: Mathematics',
      'License :: OSI Approved :: Apache Software License'
  ]
)
