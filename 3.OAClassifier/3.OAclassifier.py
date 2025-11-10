import joblib
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, confusion_matrix, roc_curve, auc, f1_score
from sklearn.preprocessing import label_binarize

# Define paths
model_dir = "./model_selection"
prediction_base_dir = "./prediction"

# Load selected features
featurelist = pd.read_csv('./feature_selection/RF_select.txt', sep=' ')
featurelist = featurelist.iloc[:, 0].tolist()

# Get all model files in the directory
model_files = [f for f in os.listdir(model_dir) if f.endswith('.pkl')]

# Define subtype mapping
subtype_mapping = {0: 'TC1', 1: 'TC2', 2: 'TC3', 3: 'TC4', 4: 'TC5', 5: 'TC6'}

# Process each model
for model_file in model_files:
    # Extract algorithm name from filename
    algorithm = model_file.replace('_best_model.pkl', '').replace('best_', '').upper()
    
    # Create directories
    prediction_alg_dir = os.path.join(prediction_base_dir, algorithm)
    os.makedirs(prediction_alg_dir, exist_ok=True)
    
    print(f"\nProcessing {algorithm} model...")
    
    # Load model
    model_path = os.path.join(model_dir, model_file)
    model = joblib.load(model_path)
    
    # Iterate over every .txt file in the prediction folder for subtype prediction 
    for filename in os.listdir(prediction_base_dir):
        if filename.endswith('.txt'):
            file_path = os.path.join(prediction_base_dir, filename)
            
            # Load external or validation dataset
            external = pd.read_csv(file_path, sep='\t', index_col=0).T
            
            # Retain only the selected features 
            X_external = external.reindex(columns=featurelist)
            X_external = X_external.fillna(0) # fill 0 if there is no corresponding features 
            
            # Predict labels
            y_external_pred = model.predict(X_external)
            y_external_pred = pd.Series(y_external_pred).map(subtype_mapping)
            
            # Predict probabilities
            y_external_pred_proba = model.predict_proba(X_external)
            
            # Create probability columns for each class
            proba_columns = {}
            for i, class_name in subtype_mapping.items():
                proba_columns[f'Probability_{class_name}'] = y_external_pred_proba[:, i]
            
            # Save results as dataFrame
            pre_results = pd.DataFrame({
                'Sample_ID': X_external.index,
                'Predicted_Label': y_external_pred,
                **proba_columns  # Unpack all probability columns
            })

            # Derive output filename from the original file name
            base_filename = os.path.splitext(filename)[0]
            output_filename = f"{base_filename}_predictions.txt"
            
            # Save the results as a tab-delimited text file
            pre_results.to_csv(
                os.path.join(prediction_alg_dir, output_filename),
                sep='\t', index=False, header=True
            )
            
            print(f"Predicted using {filename} → saved as {output_filename}")