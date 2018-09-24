import os
import sys

try:
    from setuptools import setup
except ImportError:
    sys.exit('ERROR: setuptools is required.\n')


try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements
# try:
#     from pip.req import parse_requirements
# except ImportError:
#     sys.exit('ERROR: pip is required.\n')


if os.environ.get('READTHEDOCS', None):
    # Set empty install_requires to get install to work on readthedocs
    install_requires = []
else:
    if sys.version_info[0] > 2:
        req_file = 'requirements.txt'
    else:
        req_file = 'requirements2.txt'
    try:
        reqs = parse_requirements(req_file, session=False)
    except TypeError:
        reqs = parse_requirements(req_file)
    install_requires = [str(r.req) for r in reqs]

# read version
exec(open('abutils/version.py').read())

config = {
    'description': 'Utilities for analysis of antibody NGS data',
    'author': 'Bryan Briney',
    'url': 'https://www.github.com/briney/abutils',
    'author_email': 'briney@scripps.edu',
    'version': __version__,
    'install_requires': install_requires,
    'packages': ['abutils'],
    'scripts': [],
    'name': 'abutils',
    'include_package_data': True
}

setup(**config)
