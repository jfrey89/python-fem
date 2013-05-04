try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'fempy - Finite Element Method (for) Python',
    'author': 'Will Frey',
    'url': 'URL to get it at.',
    'download_url': 'Where to download it.',
    'author_email': 'My email.',
    'version': '0.1',
    'install_requires': ['nose',
                         'numpy',
                         'scipy',
                         'matplotlib'
                         ],
    'packages': ['fempy'],
    'scripts': [],
    'name': 'python-fempy'
}

setup(**config)
