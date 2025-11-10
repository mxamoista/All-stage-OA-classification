import joblib
import pandas as pd
import numpy as np
import os
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from xgboost import XGBRegressor
from sklearn.model_selection import LeaveOneOut, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import make_scorer, accuracy_score

# Load data
X_train = pd.read_csv('./X_train.txt', sep='\t', index_col=0)
y_train = pd.read_csv('./y_train.txt', sep='\t', index_col=0)
y_train = y_train.values.ravel()  # Convert to 1D array

# Load selected features
# selected_features = pd.read_csv('./feature_selection/top_150_features.txt', sep=',', header=None,skiprows=1)[1].tolist()
selected_features = pd.read_csv('./feature_selection/RF_select.txt', sep=' ')
# selected_features = selected_features[selected_features.iloc[:, 1]>1].iloc[:, 0].tolist() #Selected features occur more than once
selected_features = selected_features.iloc[:, 0].tolist()
X_train_selected = X_train.loc[:, selected_features]

# Define models and their parameter grids
models = {
    'SVM': {
        'model': Pipeline([
            ('scaler', StandardScaler()),
            ('clf', SVC(kernel='linear', probability=True))
        ]),
        'params': {
            'clf__C': [0.1, 1, 10, 100],
            'clf__kernel': ['linear', 'rbf'],
            'clf__gamma': ['scale', 'auto']
        }
    },
    'Random Forest': {
        'model': RandomForestClassifier(random_state=42),
        'params': {
            'n_estimators': [50, 100, 200],
            'max_depth': [None, 10, 20],
            'min_samples_split': [2, 5, 10]
        }
    },
    'MLP': {
        'model': Pipeline([
            ('scaler', StandardScaler()),
            ('clf', MLPClassifier(max_iter=1000, random_state=42))
        ]),
        'params': {
            'clf__hidden_layer_sizes': [(50,), (100,), (50, 50)],
            'clf__alpha': [0.0001, 0.001, 0.01],
            'clf__learning_rate': ['constant', 'adaptive']
        }
    }
}

# Configure LOOCV and scoring
loo = LeaveOneOut()
scorer = make_scorer(accuracy_score)

# Perform grid search for each model
results = {}
for name, config in models.items():
    print(f"\n=== Tuning {name} ===")
    
    # Setup GridSearchCV with LOOCV
    grid_search = GridSearchCV(
        estimator=config['model'],
        param_grid=config['params'],
        cv=loo,
        scoring=scorer,
        n_jobs=-1,  # Use all available cores
        verbose=1
    )
    
    # Perform grid search
    grid_search.fit(X_train_selected, y_train)
    
    # Store results
    results[name] = {
        'best_model': grid_search.best_estimator_,
        'best_params': grid_search.best_params_,
        'best_score': grid_search.best_score_,
        'cv_results': grid_search.cv_results_
    }
    
    # Print summary
    print(f"Best parameters: {grid_search.best_params_}")
    print(f"Best LOOCV accuracy: {grid_search.best_score_:.4f}")
    
    # Save the best model
    output_dir = '/data2/huangbowei/radiomics/xray/SVM__RFE/model_selection'
    model_filename = os.path.join(output_dir, f'best_{name.lower().replace(" ", "_")}_model.pkl')
    joblib.dump(grid_search.best_estimator_, model_filename)


# Save all results to files
with open('./model_selection/best_parameters_summary.txt', 'w') as f:
    f.write("=== Best Parameters Summary ===\n\n")
    
    # Write each model's best parameters
    for name in results:
        f.write(f"Model: {name}\n")
        f.write(f"Accuracy: {results[name]['best_score']:.4f}\n")
        f.write("Parameters:\n")
        for param, value in results[name]['best_params'].items():
            f.write(f"  {param}: {value}\n")
        f.write("\n" + "-"*50 + "\n\n")
    

# Save all results to CSV
results_df = pd.DataFrame({
    'Model': results.keys(),
    'Best Accuracy': [results[x]['best_score'] for x in results],
    'Best Parameters': [results[x]['best_params'] for x in results]
})
results_df.to_csv('./model_selection/grid_search_results.csv', index=False)
