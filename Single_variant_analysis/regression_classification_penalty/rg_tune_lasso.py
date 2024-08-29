from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.linear_model import Lasso
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# Load the phenotype data
phenotype_name = 'merged_data_with_phenotypes.csv'  
phenotype_file = pd.read_csv(phenotype_name)
phenotype_file = phenotype_file.dropna()
# Extract phenotype (target variable)
# Selecting the first 230 columns from the merged DataFrame
selected_columns = list(phenotype_file.columns[4:230])

# Adding the specific columns requested
additional_columns = ['age', 'sex', 'BMI', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']

# Ensure we are selecting the correct columns from the merged DataFrame
final_columns = selected_columns + additional_columns

# Creating the final DataFrame with the selected columns
final_df = phenotype_file[final_columns]
y = phenotype_file['RG_ln']

# Extract genotype data (features)
X = final_df # Update as needed

# Remove rows with NaN in the target variable
mask = ~y.isna()
X = X[mask]
y = y[mask]

# Standardize the data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)
lasso = Lasso()
kf = KFold(n_splits=5, shuffle=True, random_state=42)
# Define a range of alpha values for Lasso
alpha_values = np.logspace(-10, -2, 100)
# Store results
mse_scores = []

# Loop through each alpha, perform LOOCV, and store the mean MSE score (positive value)
for alpha in alpha_values:
    lasso = Lasso(alpha=alpha, max_iter=10000, tol=0.01)
    # Negate the negative MSE scores to get positive MSE scores
    mse = -cross_val_score(lasso, X_train, y_train, cv=kf, scoring='neg_mean_squared_error').mean()
    mse_scores.append(mse)

# Find the best alpha and corresponding MSE score
best_alpha = alpha_values[np.argmin(mse_scores)]
best_mse = min(mse_scores)

# Fit the Lasso model with the best alpha on the entire dataset
best_lasso = Lasso(alpha=best_alpha, max_iter=10000, tol=0.01)
best_lasso.fit(X_train, y_train)

# Predict on the test set
y_pred = best_lasso.predict(X_test)

# Calculate and print the mean squared error for the test set
test_mse = mean_squared_error(y_test, y_pred)
r2_lasso = r2_score(y_test, y_pred)
# Store the best alpha and its score in a file
with open('rg_best_alpha_mse.txt', 'w') as f:
    f.write(f"Best Alpha: {best_alpha}\n")
    f.write(f"Best MSE Score: {best_mse}\n")
    f.write(f"Test Mean Squared Error: {test_mse}\n")
    f.write(f"R^2 score: {r2_lasso}\n")

# Plotting the MSE as a function of alpha
plt.figure(figsize=(10, 6))
plt.plot(alpha_values, mse_scores, label='Cross-Validation MSE')
plt.xscale('log')
plt.xlabel('Alpha')
plt.ylabel('MSE')
plt.title('Cross-Validation Error as a Function of Alpha')

# Highlight the minimum MSE
plt.axvline(x=best_alpha, color='r', linestyle='--', label=f'Best Alpha: {best_alpha:.2e}\nMin MSE: {best_mse:.4f}')
plt.legend()

plt.annotate(f'Best Alpha: {best_alpha:.2e}\nMSE: {best_mse:.4f}',
             xy=(best_alpha, best_mse),
             xytext=(best_alpha, best_mse + 0.5 * (max(mse_scores) - best_mse)),  # Adjust text position
             arrowprops=dict(facecolor='black', arrowstyle='->', connectionstyle='arc3'),
             horizontalalignment='right',
             verticalalignment='top')

# Save the plot to a file
plt.savefig('rg_cross_validation_mse_plot.png')

# Close the plot to free memory
plt.close()


# Get the coefficients from the best LASSO model
lasso_coefficients = best_lasso.coef_

# Create a DataFrame to hold the SNP names and their corresponding coefficients
# Create a DataFrame to hold the feature names and their corresponding coefficients
coeff_df = pd.DataFrame({'SNP': X.columns, 'Coefficient': lasso_coefficients})

# Identify non-zero coefficients
non_zero_coefficients = coeff_df[coeff_df['Coefficient'] != 0]

# Save the top 10 important variants to a file
non_zero_coefficients.to_csv('rg_non_zero.csv', index=False)
