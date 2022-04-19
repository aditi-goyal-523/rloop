import sys
import pandas as pd
import os

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

#-----------------------------------------------------------------------------------------------


#read in directory
dir=sys.argv[1]

#set output directory

save_path='processed_narrowpeak/'

#read in files
files=sorted(os.listdir(dir))
for filename in files: #for each file in dir
	print(filename)
	
	vals=filename.split('.') #split into filename, rep, pos(neg), junk, junk, narrowPeak
	root=vals[0] #name is first element from split
	#note files are processed alphabetically
	#neg will always be directly before pos
	if vals[2]=="neg": #if dealing with neg file
		negname=vals[0]+".neg.txt"
		negfile=os.path.join(dir, negname) #neg file is the file
		posname=vals[0]+".pos.txt" #temporary name for positive file
		posfile=os.path.join(dir, posname) #assign pos file to that file
	
		#read in file
		negdf=readnarrowpeak(negfile)
		posdf=readnarrowpeak(posfile)
		
		#abridge the dataframe to only contain useful cols
		negdf=abridgenarrowpeak(negdf)
		posdf=abridgenarrowpeak(posdf)
		
		#add new columns and process the data
		negdf=processnarrowpeak(negdf)
		posdf=processnarrowpeak(posdf)
		
        # #get means of the length and pileup count
		# negmeans=statsnarrowpeak(negdf)
        # posmeans=statsnarrowpeak(posdf)
	
		#write data to a new file for reference
		negout=os.path.join(save_path, negname)
		posout=os.path.join(save_path, posname)
		
		# diffname=vals[0]+".diff.txt"
		# diffout=os.path.join(save_path, diffname)
		
		with open(negout, 'w') as writer:
			writer.write(negdf.to_csv(index=False))
            #writer.write(negmeans.to_csv(index=False))
		writer.close()
		
		with open(posout, 'w') as writer:
			writer.write(posdf.to_csv(index=False))
            #writer.write(posmeans.to_csv(index=False))
		writer.close()
