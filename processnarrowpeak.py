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
    df=df[['chr', 'start', 'end', 'length', 'count', 'peak']]
    return(df)

def statsnarrowpeak(df):
    means=df[['length', 'count']].mean()
    m_length=means[0]
    m_count=means[1]

    sds=df[['length', 'count']].std()
    s_length=sds[0]
    s_count=sds[1]

    summary=pd.DataFrame({
        "length":[m_length, s_length], 
        "count":[m_count, s_count]
    })
    summary.index=["mean", "st. dev"]
    #frames = [means, sds]
    #result = pd.concat(frames)
    return(summary)

#----------------------------------------------------



df=readnarrowpeak(sys.argv[1])
df=abridgenarrowpeak(df)
df=processnarrowpeak(df)
stats=statsnarrowpeak(df)

#print out the summary statistics
print(stats)

#this section only applicable for the narrowpeak files generated from macs2 given nomenclature
#write out the new style file to a different file
sections=sys.argv[1].split(".")
df.to_csv(sections[0]+"."+sections[2]+".processed.txt", sep="\t", index=False) #includes label + pos/neg + processed.txt

	

