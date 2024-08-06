# BMN_v1.2.py
An improved python script based on a previously developed R script in biochemometrics to find potential bioactive compounds in a bioactive mixture.

See the publication for the previously developed R script: Nothias, Louis-Félix, et al. "Bioactivity-based molecular networking for the discovery of drug leads in natural product bioassay-guided fractionation." Journal of natural products 81.4 (2018): 758-767.
See the publication for the previously developed R script (hereafter referred as BMN_v1.1.r): Nothias, Louis-Félix, et al. "Bioactivity-based molecular networking for the discovery of drug leads in natural product bioassay-guided fractionation." Journal of natural products 81.4 (2018): 758-767.
Please read https://github.com/DorresteinLaboratory/Bioactive_Molecular_Networks to know the backgound of using this script.
## Usage
Example usage: python BMN_v1.2.py -f path_to_quant_activ.csv -n N -sn SN -t cor_threshold -p P_level -use_adjusted_p_value

Run "python BMN_v1.2.py -h" in the terminal to find out more about the arguments.
## Changes
Main changes compared to BMN_v1.1.r:
1) subsampling. The chemical profile we create from column chromatography would have different degrees of overlapping depending on the elution conditions. Subsampling is conducted by the user-defined parameter, N, to have better correlation analysis. Considering that the actual overlapping degree is hard to be accurately estimated, three subsampling window lengths are included (N-1, N, N+1), giving three pre-lists that are used to produce a common list.
2) data shuffling. True positive should not be influenced by the sample order. To achieve this purpose, we included another parameter, SN, to allow users to shuffle the dataframe for SN times to avoid randomness.

Other small changes:
1) user-defined correlation threshold.
2) user-defined significance level.
## Workflow
A schematic workdlow of BMN_v1.2.py:

![image](https://github.com/user-attachments/assets/4dc3347c-fb2c-4533-97b0-375394342ac7)

