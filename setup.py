from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='memeparse',
      version=version,
      description="For parsing meme files",
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='',
      author='Paul Agapow',
      author_email='paul@agapow.net',
      url='http://www.agapow.net/software/memeparse',
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
