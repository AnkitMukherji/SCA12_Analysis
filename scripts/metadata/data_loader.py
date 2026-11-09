import os
import pandas as pd
import numpy as np

def load_and_clean_data(file_path=None):
    """
    Loads the SCA12 metadata Excel sheet and performs necessary data cleaning.
    """
    if file_path is None:
        file_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'metadata', 'metadata_sca12.xlsx')
    
    # Load raw data
    df = pd.read_excel(file_path)
    
    # Clean column names (strip spaces)
    df.columns = [col.strip() for col in df.columns]
    
    # Clean Gender column (standardize capitalization and strip spaces)
    df['Gender'] = df['Gender'].astype(str).str.strip().str.capitalize()
    
    # Check for anomalous row (Record ID 60: Age=52, Disease duration=52, AO=NaN)
    # Since Age of Onset cannot be 0 for SCA12, we treat this row's AO and Disease duration
    # as missing/invalid for analyses requiring onset or duration relationships.
    # We will mark it explicitly.
    df.loc[df['Record ID'] == 60, 'AO'] = np.nan
    df.loc[df['Record ID'] == 60, 'Disease duration'] = np.nan
    
    # Clean participant ID (strip spaces)
    df['Participant ID'] = df['Participant ID'].astype(str).str.strip()
    
    return df

def get_clinical_data_for_analysis(df):
    """
    Filters rows with valid clinical scores for correlation studies.
    """
    clinical_cols = ['SARA Score', 'INAS Score', 'FTMTRS Grand Total', 'TETRAS Grand Total', 'ADL Total']
    return df.dropna(subset=clinical_cols)

def get_genetic_data_for_analysis(df):
    """
    Filters rows with valid genetic repeat sizes.
    """
    return df.dropna(subset=['LOWER ALLELE', 'UPPER ALLELE'])

if __name__ == '__main__':
    df = load_and_clean_data()
    print("Data loaded successfully!")
    print(f"Shape: {df.shape}")
    print(f"Genders: {df['Gender'].value_counts().to_dict()}")
    print(f"Anomalous Record ID 60 values after cleaning: Age={df.loc[df['Record ID'] == 60, 'Age'].values[0]}, AO={df.loc[df['Record ID'] == 60, 'AO'].values[0]}, Duration={df.loc[df['Record ID'] == 60, 'Disease duration'].values[0]}")
