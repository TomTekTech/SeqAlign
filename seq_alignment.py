import pandas as pd
import argparse as ap

parser = ap.ArgumentParser(prog='Sequence_Alignment',description='This program aligns two sequences using the needleman-wunsch algorithm',)

parser.add_argument('-f','--file', dest='sequenceFile',type=str,action='store',required=True,help='A text file contaiing the two sequences to be aligned')
parser.add_argument('-d', '--gap_penalty',dest='gapPenalty',action='store',default = -1,type = int,help='A negative integer to determine the gap penalty for the alignment (default = -1)')

args = parser.parse_args()
f = open(args.sequenceFile,'r')
a = f.readline().strip()
b = f.readline().strip()
d = args.gapPenalty

def dataframe_formation (a,b):
    seqA = a
    seqB = b

    df = pd.DataFrame()

    for i in range(0,len(seqA)+1):
        df.loc[i,0] = i * d
    
    for j in range(0,len(seqB)+1):
        df.loc[0,j] = j * d

    for i in range(1,len(seqA)+1):
        for j in range(1,len(seqB)+1):
            delete = df.loc[i-1,j] + d
            insert = df.loc[i,j-1] + d
            if seqA[i-1] == seqB[j-1]:
                m = df.loc[i-1,j-1] + 1 
            else:
                m = df.loc[i-1,j-1] - 1

            df.loc[i,j] = max(m,delete,insert)
    
    return(df)

def alignment():
    df = dataframe_formation(a,b)

    alignmentA = ''
    alignmentB = ''

    i = len(a)
    j = len(b)

    while i > 0 or j > 0:

        delete = df.iloc[i-1,j] + d
        insert = df.iloc[i,j-1] + d
        if a[i-1] == b[j-1]:
            match_score = df.iloc[i-1,j-1] + 1
        else:
            match_score = df.iloc[i-1,j-1] - 1

        if df.iloc[i,j] == match_score:
            alignmentA = a[i-1] + alignmentA
            alignmentB = b[j-1] + alignmentB
            i -= 1
            j -= 1

        elif i > 0 and df.iloc[i,j] == delete:
            alignmentA = a[i-1] + alignmentA
            alignmentB = '_' + alignmentB
            i -= 1 
        
        elif j > 0 and df.iloc[i,j] == insert:
            alignmentB = b[j-1] + alignmentB
            alignmentA = '_' + alignmentA
            j -= 1 

    print('\n')
    print(f'Sequence 1: {alignmentA}')
    print(f'Sequence 2: {alignmentB}')

alignment()

