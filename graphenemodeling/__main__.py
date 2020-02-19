#__main__.py

import os, sys

from graphenemodeling import __version__

def main():
    greeting='''
    Welcome to graphenemodeling version %s!

    For issues, contact Gregory Holdman at gholdman@protonmail.com
    or submit an issue or pull request at github.
    ''' % (__version__)

    if len(sys.argv)==1:
        print(greeting)
    if len(sys.argv)==2:
        files=os.listdir(os.path.dirname(__file__))
        if sys.argv[1]=='overview':
            from graphenemodeling import overview
            overview
        else:
            print('Unknown command: %s' % (sys.argv[1]))
        
if __name__=="__main__":
    main()
