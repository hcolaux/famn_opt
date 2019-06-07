# -*- coding: utf-8 -*-
"""
Setup and register package in PyPI
"""

from setuptools import setup

# For testing:
#
# python setup.py sdist bdist_wheel
# python -m pip install --user --upgrade twine
# python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/* -u hcolaux -p Ls2rpAn!
# twine upload dist/* -u hcolaux -p Ls2rpAn!



setup(name='famn_opt',
      version='0.3.2',
      description='Python code to optimise FAM-N pulses for MQMAS experiments',
      long_description = open('Readme.txt').read(),
      url='https://github.com/hcolaux/famn_opt',
      author='hcolaux',
      author_email='henri.colaux@normalesup.org',
      packages=['famn_opt'],
      zip_safe=False,
      license='GNU GPL',
      classifiers=[
        "Programming Language :: Python :: 3.7",
        "Development Status :: 3 - Alpha",
        "License :: GNU GPL",
        "Operating System :: OS Independent",],
      include_package_data=True,
      install_requires=[ #'sys',
                        #'os',
                        #'subprocess',
                        #'pickle',
                        'numpy',
                        'lxml',
                      'matplotlib',
                      'PyQt5',
                      ],
         extras_requires = ['pickle']
      )

