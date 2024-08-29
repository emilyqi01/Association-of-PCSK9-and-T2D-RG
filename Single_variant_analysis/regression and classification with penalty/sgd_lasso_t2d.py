from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import roc_auc_score, f1_score, make_scorer, confusion_matrix, precision_score, recall_score
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.dummy import DummyClassifier


# Load the phenotype data
phenotype_name = 'merged_data_with_phenotypes.csv'
phenotype_file = pd.read_csv(phenotype_name)
phenotype_file = phenotype_file.dropna()

# Extract phenotype (target variable)
y = phenotype_file['T2D']

# Selecting the first 230 columns from the merged DataFrame
selected_columns = list(phenotype_file.columns[4:230])

# Adding the specific columns requested
additional_columns = ['age', 'sex', 'BMI', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']

# Ensure we are selecting the correct columns from the merged DataFrame
final_columns = selected_columns + additional_columns

# Creating the final DataFrame with the selected columns
final_df = phenotype_file[final_columns]

# Extract genotype data (features)
X = final_df

# Remove rows with NaN in the target variable
mask = ~y.isna()
X = X[mask]
y = y[mask]

# Standardize the data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

# Define a narrower range of alpha values around the best one identified
alpha_values = np.logspace(-4, -1, 5)  
roc_auc_scores = []
f1_scores = []

# KFold cross-validation setup with fewer splits
kf = KFold(n_splits=3, shuffle=True, random_state=42) 

# Loop through each alpha value, perform cross-validation, and store the mean ROC AUC and F1 scores
for alpha in alpha_values:
    sgd_lasso = SGDClassifier(loss='log_loss', penalty='l1', alpha=alpha, max_iter=500, tol=1e-2, class_weight='balanced', warm_start=True)
    # Calculate ROC AUC and F1 score using cross-validation
    roc_auc = cross_val_score(sgd_lasso, X_train, y_train, cv=kf, scoring='roc_auc', n_jobs=-1).mean()
    f1 = cross_val_score(sgd_lasso, X_train, y_train, cv=kf, scoring=make_scorer(f1_score), n_jobs=-1).mean()
    roc_auc_scores.append(roc_auc)
    f1_scores.append(f1)

# Find the best alpha and corresponding ROC AUC and F1 scores
best_alpha_idx = np.argmax(roc_auc_scores)  # Optimize for ROC AUC primarily
best_alpha = alpha_values[best_alpha_idx]
best_roc_auc = roc_auc_scores[best_alpha_idx]
best_f1 = f1_scores[best_alpha_idx]

# Fit the SGDClassifier with the best alpha on the entire training set with balanced class weights
best_sgd_lasso = SGDClassifier(loss='log_loss', penalty='l1', alpha=best_alpha, max_iter=500, tol=1e-2, class_weight='balanced', warm_start=True)
best_sgd_lasso.fit(X_train, y_train)

# Predict on the test set
y_pred_proba = best_sgd_lasso.decision_function(X_test)
y_pred = best_sgd_lasso.predict(X_test)

# Calculate and print the evaluation metrics for the test set
test_roc_auc = roc_auc_score(y_test, y_pred_proba)
test_f1 = f1_score(y_test, y_pred)
test_precision = precision_score(y_test, y_pred)
test_recall = recall_score(y_test, y_pred)
conf_matrix = confusion_matrix(y_test, y_pred)

print(f'Test ROC AUC: {test_roc_auc}')
print(f'Test F1 Score: {test_f1}')
print(f'Test Precision: {test_precision}')
print(f'Test Recall: {test_recall}')
print(f'Confusion Matrix:\n{conf_matrix}')
print(f'Best alpha: {best_alpha}')
print(f'Best ROC AUC Score on Cross-Validation: {best_roc_auc}')
print(f'Best F1 Score on Cross-Validation: {best_f1}')

# Plotting the ROC AUC as a function of alpha
plt.figure(figsize=(10, 6))
plt.plot(alpha_values, roc_auc_scores, label='Cross-Validation ROC AUC')
plt.plot(alpha_values, f1_scores, label='Cross-Validation F1 Score', linestyle='--')
plt.xscale('log')
plt.xlabel('Alpha')
plt.ylabel('Score')
plt.title('Cross-Validation ROC AUC and F1 Score as a Function of Alpha')

# Highlight the best ROC AUC and F1 Score
plt.axvline(x=best_alpha, color='r', linestyle='--', label=f'Best alpha: {best_alpha:.2e}\nMax ROC AUC: {best_roc_auc:.4f}\nMax F1 Score: {best_f1:.4f}')
plt.legend()

plt.annotate(f'Best alpha: {best_alpha:.2e}\nROC AUC: {best_roc_auc:.4f}\nF1 Score: {best_f1:.4f}',
             xy=(best_alpha, best_roc_auc),
             xytext=(best_alpha, best_roc_auc - 0.05),  # Adjust text position
             arrowprops=dict(facecolor='black', arrowstyle='->', connectionstyle='arc3'),
             horizontalalignment='right',
             verticalalignment='top')

plt.show()

# Fit the SGDClassifier with the best alpha on the entire training set with balanced class weights
best_sgd_lasso = SGDClassifier(loss='log_loss', penalty='l1', alpha=best_alpha, max_iter=1000, tol=1e-3, class_weight='balanced', warm_start=True)
best_sgd_lasso.fit(X_train, y_train)

# Predict on the test set
y_pred_proba = best_sgd_lasso.decision_function(X_test)
y_pred = best_sgd_lasso.predict(X_test)

# Calculate and print the evaluation metrics for the test set
test_roc_auc = roc_auc_score(y_test, y_pred_proba)
test_f1 = f1_score(y_test, y_pred)
test_precision = precision_score(y_test, y_pred)
test_recall = recall_score(y_test, y_pred)
conf_matrix = confusion_matrix(y_test, y_pred)

print(f'Test ROC AUC: {test_roc_auc}')
print(f'Test F1 Score: {test_f1}')
print(f'Test Precision: {test_precision}')
print(f'Test Recall: {test_recall}')
print(f'Confusion Matrix:\n{conf_matrix}')
print(f'Best alpha: {best_alpha}')
print(f'Best ROC AUC Score on Cross-Validation: {best_roc_auc}')
print(f'Best F1 Score on Cross-Validation: {best_f1}')

coefficients = best_sgd_lasso.coef_.flatten()

# Extract the non-zero coefficients
non_zero_coefficients = coefficients[coefficients != 0]
non_zero_features = np.array(final_columns)[coefficients != 0]

# Create a DataFrame to display non-zero coefficients with corresponding features
non_zero_df = pd.DataFrame({
    'Feature': non_zero_features,
    'Coefficient': non_zero_coefficients
})

# Sort the DataFrame by the absolute value of the coefficients in descending order
non_zero_df = non_zero_df.reindex(non_zero_df.Coefficient.abs().sort_values(ascending=False).index)

print("\nSignificant Features (Non-Zero Coefficients):")
print(non_zero_df)

# Create a DummyClassifier as a baseline
dummy_clf = DummyClassifier(strategy='most_frequent')
dummy_clf.fit(X_train, y_train)

# Predict on the test set using the DummyClassifier
y_dummy_pred_proba = dummy_clf.predict_proba(X_test)[:, 1]

# Calculate the baseline ROC AUC
baseline_roc_auc = roc_auc_score(y_test, y_dummy_pred_proba)
print(f'Baseline ROC AUC (Dummy Classifier): {baseline_roc_auc}')

# Fit the SGDClassifier with the best alpha on the entire training set with balanced class weights
best_sgd_lasso = SGDClassifier(loss='log_loss', penalty='l1', alpha=best_alpha, max_iter=2000, tol=1e-3, class_weight='balanced', warm_start=True)
best_sgd_lasso.fit(X_train, y_train)

# Predict on the test set
y_pred_proba = best_sgd_lasso.decision_function(X_test)
y_pred = best_sgd_lasso.predict(X_test)

# Calculate and print the evaluation metrics for the test set
test_roc_auc = roc_auc_score(y_test, y_pred_proba)
print(f'Test ROC AUC (SGD Classifier): {test_roc_auc}')

# Compare the baseline ROC AUC with the model's ROC AUC
improvement = test_roc_auc - baseline_roc_auc
print(f'Improvement in ROC AUC over baseline: {improvement}')

# (Continue with the extraction of non-zero coefficients if needed)

from sklearn.dummy import DummyRegressor
from sklearn.metrics import mean_squared_error, r2_score
# Create a DummyRegressor as a baseline
dummy_regressor = DummyRegressor(strategy='mean')
dummy_regressor.fit(X_train, y_train)

# Predict on the test set using the DummyRegressor
y_dummy_pred = dummy_regressor.predict(X_test)

# Calculate the baseline MSE
baseline_mse = mean_squared_error(y_test, y_dummy_pred)
print(f'Baseline MSE (Dummy Regressor): {baseline_mse}')

r2_dummy = r2_score(y_test, y_dummy_pred)
print(f'Baseline R^2 (Dummy Regressor): {r2_dummy}')