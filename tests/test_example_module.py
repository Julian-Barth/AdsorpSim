import pytest
import pandas as pd
import os  # Import os module to resolve NameError

# Function to read the CSV file, ensure it exists before reading
def download_data(file_name):
    file_path = file_name
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    else:
        # Reading the CSV file with proper separator
        return pd.read_csv(file_path, sep=";")

# Test to check if the file exists and is read correctly
def test_download_data():
    file_path = "data/Dataset PPC.csv"
    
    # Check if the file exists
    if not os.path.exists(file_path):
        pytest.fail(f"Test file '{file_path}' not found.")
    
    df = download_data(file_path)
    assert not df.empty, "Dataframe is empty after reading the file."
    assert 'Adsorbent' in df.columns, "'Adsorbent' column is missing from the data."

# Test to check if adding an adsorbent works
def test_add_adsorbent():
    # Load the CSV file into a dataframe
    df = download_data("data/Dataset PPC.csv")
    
    # Check initial number of adsorbents
    initial_count = len(df)
    
    # Add a new adsorbent
    new_adsorbent = "NewAdsorbent"
    new_entry = {"Adsorbent": new_adsorbent, "OtherColumn": "SomeValue"}  # Adjust as needed
    df = df.append(new_entry, ignore_index=True)
    
    # Check if the new adsorbent has been added
    assert len(df) == initial_count + 1, "Adsorbent was not added."
    assert new_adsorbent in df['Adsorbent'].values, "New adsorbent is not found in the dataframe."

# Test to check if adding an existing adsorbent raises an error
def test_add_existing_adsorbent():
    df = download_data("data/Dataset PPC.csv")
    
    # Try adding an already existing adsorbent (assuming "ExistingAdsorbent" is in the CSV file)
    with pytest.raises(ValueError):
        df = df.append({"Adsorbent": "ExistingAdsorbent", "OtherColumn": "SomeValue"}, ignore_index=True)

# Test case to ensure the data is loaded properly and contains expected columns
def test_data_columns():
    df = download_data("data/Dataset PPC.csv")
    expected_columns = ["Adsorbent", "OtherColumn"]  # Adjust this as per your actual columns
    for column in expected_columns:
        assert column in df.columns, f"'{column}' column is missing in the data."
