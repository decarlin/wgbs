import numpy as np
import scanpy as sc
import os
import pandas as pd
import pybedtools
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression


metadata=pd.read_csv('/data/combined_wgbs/combined_metadata.csv', index_col=0)

os.chdir('/data/combined_wgbs/')
mch_adata = sc.read_h5ad(filename='/data/combined_wgbs/adrc_chrom1k_CGN_01var.h5ad')

#pc2_loadings = mch_adata.obsm['X_pca'][:, 1]                                                                                                                
#max_index = np.argmax(pc2_loadings)                                                                                                                         
#print("Max PC1 loading:", pc1_loadings[max_index])                                                                                                          
#print("Corresponding sample:", adata.obs_names[max_index])                                                                                                  

#pc2_series = pd.Series(pc2_loadings, index=mch_adata.obs_names)                                                                                             

# Sort descending to get top                                                                                                                                 
#n=4                                                                                                                                                         
#top_pc2 = pc2_series.sort_values(ascending=False).head(n)                                                                                                   

#blacklist with the cancer samples                                                                                                                           

blacklist = ['ADRC3440-16','PDBP16574','PDBP16676','PDBP16551','ADRC3533-16',
             'PDBP1242', 'PDBPDLS067S', 'PDBPDLS072S', 'PDBPDLS362S', 'PDBPDLS403S',
       'PDBPDLS459S', 'PDBPDLS631S', 'PDBPDLS720S', 'PDBPDLS997S',
       'ADRCDLS072P', 'ADRCDLS459P', 'ADRCDLS720P', 'ADRCRep393P']
mch_adata = mch_adata[~mch_adata.obs_names.isin(blacklist), :]

mch_adata.obs.index = mch_adata.obs.index.str.replace('_', '-')
mch_adata.obs=mch_adata.obs.join(metadata)

nan_indices = mch_adata.obs[pd.isna(mch_adata.obs['sex'])].index

#start with pdbp

X = mch_adata.X[mch_adata.obs['study']=='ADRC', :]
DX_labels = mch_adata.obs["dx"][mch_adata.obs['study']=='ADRC'].values

# Identify unique DX categories (excluding "CTRL")                                                                                                           
categories = [dx for dx in np.unique(DX_labels) if (dx != "CTRL" and dx != "Unknown")]

# Set up 5-fold cross-validation                                                                                                                             
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

plt.figure(figsize=(8, 6))  # Initialize ROC plot

# Iterate through each category vs. "CTRL"

# Number of top features to keep
K = 500

for category in categories:
    print(f"Processing {category} vs. CTRL...")
    # Binary classification: 1 = category, 0 = "CTRL"                                                                                                        
    y = np.array([1 if dx == category else 0 if dx == "CTRL" else -1 for dx in DX_labels])
    mask = y != -1
    X_filtered, y_filtered = X[mask], y[mask]
    # Storage for ROC data                                                                                                                                   
    mean_fpr = np.linspace(0, 1, 100)  # Standardized FPR points                                                                                             
    tprs, aucs = [], []
    # Perform 5-fold CV                                                                                                                                      
    for train_idx, test_idx in cv.split(X_filtered, y_filtered):
        X_train, X_test = X_filtered[train_idx], X_filtered[test_idx]
        y_train, y_test = y_filtered[train_idx], y_filtered[test_idx]
        # Feature selection: t-test (ANOVA F-test) on training data only
        #selector = SelectKBest(score_func=f_classif, k=K)
        #X_train_sel = selector.fit_transform(X_train, y_train)
        #X_test_sel = selector.transform(X_test)  # Apply same selection to test
        # Train model
        model = RandomForestClassifier(n_estimators=50)
        model.fit(X_train, y_train)
        # Predict probabilities for class 1
        y_probs = model.predict_proba(X_test)[:, 1]
        # Compute ROC curve and AUC
        fpr, tpr, _ = roc_curve(y_test, y_probs)
        roc_auc = auc(fpr, tpr)
        # Interpolate TPR for averaging
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(roc_auc)
    # Compute mean and std for TPRs                                                                                                                          
    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)
    # Compute mean AUC                                                                                                                                       
    mean_auc = np.mean(aucs)
    std_auc = np.std(aucs)
    # Plot mean ROC curve with confidence band                                                                                                               
    plt.plot(mean_fpr, mean_tpr, label=f"{category} (AUC = {mean_auc:.2f} ± {std_auc:.2f})")
    plt.fill_between(mean_fpr, mean_tpr - std_tpr, mean_tpr + std_tpr, alpha=0.2)

# Plot chance level                                                                                                                                          
plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Chance")

# Final plot settings                                                                                                                                        
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Cross-Validated ROC Curves for DX Categories vs. CTRL")
plt.legend()
plt.show()

# Plot chance level                                                                                                                                          
plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Chance")

# Final plot settings                                                                                                                                        
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Cross-Validated ROC Curves for DX Categories vs. CTRL")
plt.legend()
plt.show()

output_path = "/data/combined_wgbs/roc_curves_RF50_NoSelection_adrc.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")

#try regression

                                                                                                                                                              
del mean_auc                                                                                                                                                  
del std_auc                                                                                                                                                   
plt.clf()                                                                                                                                                     
                                                                                                                                                              
# Iterate through each category vs. "CTRL"                                                                                                                    
for category in ['DLB']:                                                                                                                                      
    print(f"Processing {category} vs. CTRL...")                                                                                                               
    # Binary classification: 1 = category, 0 = "CTRL"                                                                                                         
    y = np.array([1 if dx == category else 0 if dx == "CTRL" else -1 for dx in DX_labels])                                                                    
    mask = y != -1                                                                                                                                            
    X_filtered, y_filtered = X[mask], y[mask]                                                                                                                 
    # Storage for ROC data                                                                                                                                    
    mean_fpr = np.linspace(0, 1, 100)  # Standardized FPR points                                                                                              
    tprs, aucs = [], []                                                                                                                                       
    # Perform 5-fold CV                                                                                                                                       
    for train_idx, test_idx in cv.split(X_filtered, y_filtered):                                                                                              
        model = LogisticRegression(penalty='l1',C=2, solver='liblinear', random_state=42)                                                                     
        model.fit(X_filtered[train_idx], y_filtered[train_idx])                                                                                               
        y_probs = model.predict_proba(X_filtered[test_idx])[:, 1]  # Probabilities for class 1                                                                
        # Compute ROC curve                                                                                                                                   
        fpr, tpr, _ = roc_curve(y_filtered[test_idx], y_probs)                                                                                                
        roc_auc = auc(fpr, tpr)                                                                                                                               
        # Interpolate TPR values at common FPR points                                                                                                         
        interp_tpr = np.interp(mean_fpr, fpr, tpr)                                                                                                            
        interp_tpr[0] = 0.0                                                                                                                                   
        # Ensure the curve starts at (0,0)                                                                                                                    
        tprs.append(interp_tpr)                                                                                                                               
        aucs.append(roc_auc)                                                                                                                                  
    # Compute mean and std for TPRs                                                                                                                           
    mean_tpr = np.mean(tprs, axis=0)                                                                                                                          
    std_tpr = np.std(tprs, axis=0)                                                                                                                            
    # Compute mean AUC                                                                                                                                        
    mean_auc = np.mean(aucs)                                                                                                                                  
    std_auc = np.std(aucs)                                                                                                                                    
    # Plot mean ROC curve with confidence band                                                                                                                
    plt.plot(mean_fpr, mean_tpr, label=f"{category} (AUC = {mean_auc:.2f} ± {std_auc:.2f})")                                                                  
    plt.fill_between(mean_fpr, mean_tpr - std_tpr, mean_tpr + std_tpr, alpha=0.2)                                                                             
                                                                                                                                                              
# Plot chance level                                                                                                                                           
plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Chance")                                                                                        
                                                                                                                                                              
# Final plot settings                                                                                                                                         
plt.xlabel("False Positive Rate")                                                                                                                             
plt.ylabel("True Positive Rate")                                                                                                                              
plt.title("Cross-Validated ROC Curves for DX Categories vs. CTRL")                                                                                            
plt.legend()                                                                                                                                                  
plt.show()                                                                                                                                                    
                                                                                                                                                              
output_path = "/data/combined_wgbs/pdbp_roc_curves_L1c2logistic_pdbpDLB_noSelection.tiff"
plt.savefig(output_path, format="tiff", dpi=300, bbox_inches="tight")                                                                                         


