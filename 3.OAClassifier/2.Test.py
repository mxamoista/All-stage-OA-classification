import joblib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.metrics import accuracy_score, confusion_matrix, roc_curve, auc, f1_score
from sklearn.preprocessing import label_binarize

# Define subtype mapping
subtype_mapping = {'TC1': 0, 'TC2': 1, 'TC3': 2, 'TC4': 3, 'TC5': 4, 'TC6': 5}

# Define paths
model_dir =  "./model_selection"
picture_base_dir = "./picture"
prediction_base_dir = "./prediction"

# Load selected features
selected_features = pd.read_csv('./feature_selection/RF_select.txt', sep=' ')
# selected_features = selected_features[selected_features.iloc[:, 1]>1].iloc[:, 0].tolist() #Selected features occur more than once
selected_features = selected_features.iloc[:, 0].tolist()

# Load test data
X_test = pd.read_csv('./X_test.txt', sep='\t', index_col=0)
y_test = pd.read_csv('./y_test.txt', sep='\t', index_col=0)
X_test_selected = X_test.loc[:, selected_features]

# Get all model files in the directory
model_files = [f for f in os.listdir(model_dir) if f.endswith('.pkl')]

# Process each model
for model_file in model_files:
    # Extract algorithm name from filename
    algorithm = model_file.replace('_best_model.pkl', '').upper()
    
    # Create directories
    picture_alg_dir = os.path.join(picture_base_dir, algorithm)
    prediction_alg_dir = os.path.join(prediction_base_dir, algorithm)
    os.makedirs(picture_alg_dir, exist_ok=True)
    os.makedirs(prediction_alg_dir, exist_ok=True)
    
    print(f"\nProcessing {algorithm} model...")
    
    # Load model
    model_path = os.path.join(model_dir, model_file)
    model = joblib.load(model_path)
    
    # Predict on the test set
    y_pred = model.predict(X_test_selected)
    
    # Calculate accuracy
    accuracy = accuracy_score(y_test, y_pred)
    print(f"Accuracy on test set: {accuracy:.4f}")
    
    # Binarize the output for ROC curve
    y_test_bin = label_binarize(y_test, classes=[0, 1, 2, 3, 4, 5])
    n_classes = y_test_bin.shape[1]
    
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test_bin[:, i], model.predict_proba(X_test_selected)[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])
    
    # Save AUC values to algorithm-specific folder
    auc_file = os.path.join(prediction_alg_dir, 'auc_values.txt')
    with open(auc_file, 'w') as f:
        for i in range(n_classes):
            f.write(f"TC{i+1}: {roc_auc[i]:.4f}\n")
    
    # Calculate F1 scores for each class
    f1_scores = f1_score(y_test, y_pred, average=None)
    f1_file = os.path.join(prediction_alg_dir, 'f1_scores.txt')
    with open(f1_file, 'w') as f:
        for i, score in enumerate(f1_scores):
            f.write(f"TC{i+1}: {score:.4f}\n")
    
    # Plot individual ROC curves for each TC subtype
    colors = ['#00CED1', '#984EA3', '#FF0000', '#0000CD', '#EE82EE', '#87CEEB']
    for i, color in zip(range(n_classes), colors):
        plt.figure(figsize=(8, 8))
        plt.plot(fpr[i], tpr[i], color=color, lw=2,
                 label=f'ROC curve of TC{i+1} (area = {roc_auc[i]:0.2f})')
        plt.plot([0, 1], [0, 1], 'k--', lw=2)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('1-specificity', fontsize=12)
        plt.ylabel('Sensitivity', fontsize=12)
        plt.title(f'ROC Curve for TC{i+1}', fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.legend(loc="lower right")
        plt.savefig(os.path.join(picture_alg_dir, f'TC{i+1}_roc_curve.png'), 
                   format='png', dpi=300)
        plt.close()
    
    # Confusion Matrix with proportions (0-1.0)
    cm = confusion_matrix(y_test, y_pred)
    cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    
    plt.figure(figsize=(12, 10))
    im = plt.imshow(cm_normalized, interpolation='nearest', cmap=plt.cm.YlGnBu, vmin=0, vmax=1)
    plt.grid(False)
    cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
    cbar.set_label('Proportion', rotation=270, labelpad=20, fontsize=14)
    plt.title('Normalized Confusion Matrix', fontsize=18, pad=20)
    plt.xlabel('Predicted Label', fontsize=16)
    plt.ylabel('True Label', fontsize=16)
    
    tick_marks = np.arange(len(subtype_mapping))
    plt.xticks(tick_marks, ['TC1', 'TC2', 'TC3', 'TC4', 'TC5', 'TC6'], 
               rotation=45, fontsize=14)
    plt.yticks(tick_marks, ['TC1', 'TC2', 'TC3', 'TC4', 'TC5', 'TC6'], 
               fontsize=14)
    
    plt.tight_layout()
    plt.savefig(os.path.join(picture_alg_dir, 'confusion_matrix.png'), 
                format='png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Results saved in:\n- {picture_alg_dir}\n- {prediction_alg_dir}")

print("\nAll models processed successfully!")