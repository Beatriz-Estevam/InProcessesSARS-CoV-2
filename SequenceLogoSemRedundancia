import argparse
import os
import re
import pandas as pd
import logomaker
import numpy as np
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("BordasPath", help="arquivos com as bordas")
parser.add_argument("SequenceLogoInputFile",
                    help="arquivo a ser usado como input para o desenvolvimento do sequence logo")
parser.add_argument('SeqType', help='aa ou nt')
parser.add_argument('NumBorder', type=int, help='Número de 1-15')
args = parser.parse_args()

bordasFiles = os.scandir(args.BordasPath)
seqdict = {}
lentype = {'nt': 24,
           'aa': 8
           }
alphatype ={'nt': ['A', 'C', 'G', 'T'],
            'aa': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
            }

patternString = '.+\.bordas\.'+str(args.NumBorder)+'_\d+_\d+\.'+args.SeqType+'\.fa$'
pattern = re.compile(patternString)
for file in bordasFiles:
    match = re.search(pattern, file.name)
    if match:
        filepath = args.BordasPath+'/'+file.name
        with open(filepath, "r") as file_object:
            for line in file_object:
                line = line.rstrip()
                match2 = re.search(r'^>', line)
                if not match2:
                    if len(line.upper()) == lentype[args.SeqType]:
                        seqdict[line.upper()] = 1
                    else:
                        print('O tamanho da sequência não é o esperado', file.name)

        print(file.name)
posdict ={}
for seq in seqdict.keys():
    print(seq)
    for e in range(len(seq)):
        if e not in posdict.keys():
            posdict[e] = {}
            if seq[e] not in posdict[e]:
                posdict[e][seq[e]] = 1
        else:
            if seq[e] not in posdict[e]:
                posdict[e][seq[e]] = 1
            else:
                posdict[e][seq[e]] = posdict[e][seq[e]]+1

matrix = {'A': [],
          'C': [],
          'G': [],
          'T': []
          }
with open(args.SequenceLogoInputFile, 'w') as outputfile:
    outputfile.write(f'  \t    A\t    C\t    G\t    T\n')
    for pos in posdict.keys():
        outputfile.write(f'{pos:<2}\t')
        for res in alphatype[args.SeqType]:
            if res in posdict[pos]:
                freq = posdict[pos][res] / len(seqdict)
                matrix[res].append(freq)
                outputfile.write(f'{freq:.3f}\t')
            else:
                matrix[res].append(0)
                outputfile.write(f'0.000\t')
        outputfile.write('\n')

#SEQUENCE LOGO
DataFrame = pd.DataFrame(matrix)
fonts = logomaker.list_font_names()
fonts[0] = 'Arial Rounded MT Bold[1]'
out_df = logomaker.validate_matrix(df=DataFrame, matrix_type='probability', )
#ate aqui funcionando perfeitamente
ss_logo = logomaker.Logo(out_df,
                         width=.8,
                         vpad=.05,
                         fade_probabilities=True,
                         stack_order='small_on_top',
                         color_scheme='dodgerblue',
                         font_name=('Arial Rounded MT Bold[1]'))

ss_logo.style_spines(spines=('left', 'right'), visible=False) #color='k' é preto
ss_logo.ax.set_xticks(range(len(out_df)))
ss_logo.ax.set_xticklabels('%+d'%x for x in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                             21, 22, 23, 24]) #criar essa lista pra nt e aa
ss_logo.ax.set_yticks([0, 0.5, 1])
ss_logo.ax.set_ylabel('probability')
