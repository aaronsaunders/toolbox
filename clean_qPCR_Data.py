import fileinput
import sys

def cleanLine(line):
    elements = [ element.strip() for element in line.strip().split('\t') ]
    return '\t'.join(elements)

def cleanFname(fname):
    for line in fileinput.input(fname, inplace= True):
        print cleanLine(line)

def main(args):
    if args:
        for fname in args:
            cleanFname(fname)

    else:
        with sys.stdin as fh:
            newline = cleanLine(line)
            print newline

if __name__ == '__main__':
    main(sys.argv[1:])
