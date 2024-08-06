import pandas as pd
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr
from statsmodels.stats import multitest
import sys
import argparse

"""
This script is used to do Pearson correlation in the bioactive molecular network workflow. Originally developed by Nothias, Louis-Félix, et al. 
"Bioactivity-based molecular networking for the discovery of drug leads in natural product bioassay-guided fractionation." Journal of natural products 81.4 (2018): 758-767."

Modified by Jianzhong Zhu to allow adaption of the script to different kinds of fraction data and more accurate predictions based on common ions in different lists.

Note: with this script, one doesn't need to reset the index during data processing in MZmine.

Usage: python BMN_v1.2.py -f path_to_quant_activ.csv -n N -sn SN -t cor_threshold -p P_level -use_adjusted_p_value

E: hainu_carl88@163.com

Options:

  -h, --help               Show this help message and exit
  -v, --verbose            Enable verbose output.
  -f, --filename           Path to the filename. This is the feature_quant.csv file one can get after processing raw MS data with MZmine.      
  -n, --number             The number of overlapping features.
  -sn, --shuffling number  The number of shuffling times.
  -t, --threshold          A user-defined threshold for the correlation coefficient to narrow down the features.
  -p, --p level			   Significance level.
  -use_adjusted_p_value    A boolean flag. If present, use_adjusted_p_value will be True, and if not present, it will be False.
"""

def main():
	parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	
	# Add other command-line arguments as needed
	parser.add_argument('-f', '--filename', action='store', type=str, help='Path to the filename. This is the feature_quant.csv file one can get after processing raw MS data with MZmine.')
	parser.add_argument('-n', '--number', action='store', type=int, help='The number of overlapping features.')
	parser.add_argument('-sn', '--shuffle_number', action='store', type=int, help='The times of shuffling the samples. When n equals the number of samples, shuffling or not does not matter. But when n is lower than the number of samples, it is better to shuffle the samples to reduce randomness and increase credibility.')
	parser.add_argument('-t', '--threshold', action='store', type=float, help='A user-defined threshold for the correlation coefficient to narrow down the features.')
	parser.add_argument('-p', '--p_level', action='store', type=float, help='significance level')
	parser.add_argument('-use_adjusted_p_value', action='store_true', help='A boolean flag. If present, use_adjusted_p_value will be True, and if not present, it will be False.')
	
	args = parser.parse_args()
	
	# Check if there are enough arguments
	if len(sys.argv) < 7:
		print(f"Usage: python BMN_v1.2.py -f path_to_example_quant_activ.csv -n N -sn SN -t cor_threshold -p P_level -use_adjusted_p_value")
		sys.exit(1)
	
	return args

if __name__ == "__main__":
    args = main()
	
path = args.filename

# read file and transpose the matrix
input_tab = pd.read_csv(path)
tab1 = input_tab.iloc[1:,3:].T

# rename column names
new_column_names = input_tab['row m/z'].astype(str)+'_'+input_tab['row retention time'].astype(str)
new_column_names = new_column_names.iloc[1:]
tab1.rename(columns=dict(zip(tab1.columns, new_column_names)), inplace=True)

# insert two new columns
tab1["Sample_Names"] = tab1.index
tab1.insert(0,'Sample_Names',tab1.pop('Sample_Names'))
tab1.insert(1,'Bioactivity',input_tab.iloc[0,3:].T)

# remove row names
tab1.reset_index(drop=True,inplace=True)

# add 1 to all the peak areas
# normalization of the peak area on the rows
def custom_normalization(row):
    return (row+1)/sum(row+1)

selected_columns = tab1.columns[2:]
tab1[selected_columns] = tab1[selected_columns].apply(custom_normalization, axis=1).apply(pd.Series)

# standardize on the columns
scaler = StandardScaler()
tab1[selected_columns] = tab1[selected_columns].apply(lambda x: scaler.fit_transform(x.values.reshape(-1, 1)).flatten()) 

# define a function to do correlation test with designated number of overlapping features
def find_correlation(dataFrame, col1, col2, N):
    """
    dataFrame: a pandas dataFrame. samples on the column direction and features on the rwo direction.
    col1: string, a column of data used to do Pearson correlation test; this is the activity data.
    col2: string, a column of data used to do Pearson correlation test; this can be all the features.
    N: integer, the number of overlapping features.
    """
    
    num = len(dataFrame) - N + 1
    cor_coeff_list, p_value_list = [], []
    for i in range(num):
        cor_coeff_temp, p_value_temp = pearsonr(dataFrame[col1][i:i+N],dataFrame[col2][i:i+N])
        cor_coeff_list.append(cor_coeff_temp)
        p_value_list.append(p_value_temp)
    cor_coeff = max(cor_coeff_list)
    index = cor_coeff_list.index(cor_coeff)
    p_value = p_value_list[index]
    
    return cor_coeff, p_value, index

# define a function to get the prediction list. --> find correlation on all the features
def get_prediction_list(dataFrame, Bioactivity, N, p_value_column):
    """
    dataFrame: a pandas dataFrame. samples on the column direction and features on the rwo direction.
    Bioactivity: string, a column of data used to do Pearson correlation test; this is the activity data.
    col2: string, a column of data used to do Pearson correlation test; this can be all the features
    N: integer, the number of overlapping features.
    p_value_column: argument to use Bonferroni correction or not
    """

    # do correlation for all the features
    cor_coeffs, p_values, indexes = zip(*[find_correlation(dataFrame, Bioactivity, feature, N) for feature in dataFrame.columns[2:]])

    # attach the cor_coeffs, and p_values into the dataframe
    cor_tab = pd.DataFrame(("cor_coeffs",0)+(cor_coeffs)).T
    cor_tab.columns = tab1.columns
    p_tab = pd.DataFrame(("p_values",0)+(p_values)).T
    p_tab.columns = tab1.columns
    tab2 = pd.concat([cor_tab,p_tab,tab1])

    # transpose the table
    tab3 = tab2.T
    # use the first row as column names
    tab3.columns = tab3.iloc[0]
    # remove the first row
    tab3 = tab3.iloc[1:]
    # use the rownames as a new column for feature ID
    tab3.insert(0,"IDs",tab3.index)
    # remove the rownames so that the new index can match the index from input_tab in order to put the row ID in input_tab to tab3
    tab3 = tab3.reset_index(drop=True)
    tab3.insert(0,"row_IDs",input_tab['row ID'])
    tab3.iloc[0,0] = ''

    # apply bonferroni correction
    tab4 = tab3.loc[1:]
    tab4.insert(4,'adjusted_p_values',multitest.multipletests(tab4['p_values'], method='bonferroni')[1])

    tab4_sorted = tab4.sort_values(by=p_value_column,ascending=True)

    # putting the first row in tab3 back
    first_row = pd.DataFrame(tab3.loc[0]).T
    first_row.insert(4,'adjusted_p_value',0)
    first_row = first_row.reindex(columns=tab4_sorted.columns)

    predList = pd.concat([first_row, tab4_sorted], axis=0, ignore_index=True)

    return predList

# define a function to find the common ions in a list containing many lists
def find_common_ions(list_of_lists, cor_threshold, p_value_column):
    """
    list_of_lists: a single list of the prediction lists.
    """
    if len(list_of_lists) < 2:
        print("Error: At least two lists are required.")
        return None
    
    # Initialize common_list with the first list in the input
    common_list = list_of_lists[0]
    common_list = common_list[(common_list["cor_coeffs"] > cor_threshold) & (common_list[p_value_column] < P_level)]
    
    for x in range(1, len(list_of_lists)):
        current_list = list_of_lists[x]
        current_list = current_list[(current_list["cor_coeffs"] > cor_threshold) & (current_list[p_value_column] < P_level)]
        
        # Merge with the common_list on the "IDs" column
        common_list = pd.merge(common_list, current_list, on="IDs", how="inner")
    
    return common_list

N = int(args.number)
SN = int(args.shuffle_number)
P_level = float(args.p_level)

# decide to find common ions or not based on the size of N
p_value_column = "adjusted_p_values" if args.use_adjusted_p_value else "p_values"

# set the columns to keep to reduce redundant columns in the final list.
columns_to_keep = ["row_IDs", "IDs","cor_coeffs","p_values","adjusted_p_values"]

if N <= (tab1.shape[0] - 1):
    semifinal_lists = []
    if SN == 0:
        predList_1 = get_prediction_list(tab1, "Bioactivity", N-1, p_value_column).loc[:,columns_to_keep]
        predList_2 = get_prediction_list(tab1, "Bioactivity", N, p_value_column).loc[:,columns_to_keep]
        predList_3 = get_prediction_list(tab1, "Bioactivity", N+1, p_value_column).loc[:,columns_to_keep]
        
        all_lists = [predList_1, predList_2, predList_3]
        final_list = find_common_ions(all_lists, args.threshold, p_value_column)
		
    else:
        for i in range(SN):
            shuffled_tab1 = tab1.sample(frac=1)
            predList_1 = get_prediction_list(shuffled_tab1, "Bioactivity", N-1, p_value_column).loc[:,columns_to_keep]
            predList_2 = get_prediction_list(shuffled_tab1, "Bioactivity", N, p_value_column).loc[:,columns_to_keep]
            predList_3 = get_prediction_list(shuffled_tab1, "Bioactivity", N+1, p_value_column).loc[:,columns_to_keep]

            all_lists = [predList_1, predList_2, predList_3]
            temp_list = find_common_ions(all_lists, args.threshold, p_value_column)
            semifinal_lists.append(temp_list)

        final_list = find_common_ions(semifinal_lists, args.threshold, p_value_column) if SN > 1 else temp_list
        
else:
    final_list = get_prediction_list(tab1, "Bioactivity", N, p_value_column)
    final_list = final_list[(final_list["cor_coeffs"] > args.threshold) & (final_list[p_value_column] < P_level)]

# write the final list
final_list.to_csv(path[:-4]+'_v1.2_N'+str(N)+'_SN'+str(SN)+'_t'+str(args.threshold)+'_p'+str(args.p_level)+'.csv', index=False)

print("Success! Target file has been generated. Thanks for using this script.")