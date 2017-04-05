import pandas as pd 

def find_common_between_data_frames (df1, df2) :
	common = pd.merge(df1, df2, on='pos', how='inner')			
	return common