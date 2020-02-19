#__main__.py

import sys

from graphenemodeling import __version__

def main():
    greeting='''
    Welcome to graphenemodeling version %s!

    For issues, contact Gregory Holdman at gholdman@protonmail.com
    or submit an issue or pull request at github.
    ''' % (__version__)

    print(greeting)
if __name__=="__main__":
    main()
