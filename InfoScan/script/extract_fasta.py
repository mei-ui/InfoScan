#vim fasta-extract.py
# -*- coding: utf-8 -*-
import sys 

def usage():
    print('Usage: python script.py [fasta_file] [idlist_file] [outfile_name]')


def main():
    outf = open(sys.argv[3],'w')
    dict = {}
    with open(sys.argv[1], 'r') as fastaf:
        for line in fastaf:
            if line.startswith('>'):
                name = line.strip().split()[0][1:]
                dict[name] = ''
            else:
                dict[name] += line.replace('\n','')
     
    with open(sys.argv[2],'r') as listf:
        for row in listf:
            row = row.strip()
            for key in dict.keys():
                if key == row:
                    outf.write('>' + key+ '\n')
                    outf.write(dict[key] + '\n')
    outf.close()


try:
    main()
except IndexError:
    usage()
