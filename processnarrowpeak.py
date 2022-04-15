import sys
import pandas as pd

def readnarrowpeak(file):
    csv=pd.read_csv(file, sep="\s+", header=None)
    df=pd.DataFrame(csv)
    df.columns=['chr', 'start', 'end', 'name', 'count', 'X1', "X2", "X3", "X4", "peak"]
    return(df)

def abridgenarrowpeak(df):
    df=df[['chr', 'start', 'end', 'count', 'peak']]
    return(df)

def processnarrowpeak(df):
    df['length']=df['end']-df['start']
    df['count']=df['count']/10
    df['peak']=df['start']+df['peak']
    df=df[['start', 'end', 'length', 'count', 'peak']]
    return(df)

def statsnarrowpeak(df):
    means=df[['length', 'count']].mean()
    return(means)


#----------------------------------------------------

df=readnarrowpeak(sys.argv[1])
df=abridgenarrowpeak(df)
df=processnarrowpeak(df)

print(df)

print(statsnarrowpeak(df))