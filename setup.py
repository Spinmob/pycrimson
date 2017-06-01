

version = '0.1'



from distutils.core import setup
setup(name           = 'PyCrimson',
      version        = version,
      description    = 'Python interface for Per Vices Crimson',
      author         = 'Jack Sankey',
      author_email   = 'jack.sankey@gmail.com',
      url            = 'https://github.com/Spinmob/pycrimson/wiki',
      packages       = ['pycrimson'],
      package_dir    = {'pycrimson' : '.'})
