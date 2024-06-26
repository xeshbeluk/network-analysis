# -*- coding: utf-8 -*-
"""NE-protein mutations

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1FYwcqPIHjKs8FMTIEjmMXWvWQV5zEfWf

# **Network Exploration of Maize Protein Mutations**

## **Abstract**

The purpose of the study is to explore the relationship between the maize protein protein interaction network and the fitness cost a genetic mutation. Specifically, this study explored the betweeness centrality, degree, and average neighbor degree of specific protein coding genes targeted in an unreleased prior dataset and study.

## **Methods**

Python was the programing language used for this analysis. The igraph package was used for network characteristics exploration while statsmodels.api and scikit-learn were used for the statistical analysis.

Specifically, betweeness centrality was calulated for each of the 123 target nodes using the igraph betweeness method. Node degree was found using the igraph degree method. Lastly, average neighbor degree was found by using a custom loop that used the igraph degree method for all node neighbord and calulated the average.

Multiple regression was run using the statsmodels.api OLS method to find the p value for each coefficent and R^2 of the final model. The statsmodels.api OLS was them used to explore polynomial regression for a second, third, and fourth degree model with associated AIC and BIC measurments used for model selection.

Lastly, the model fit was explored on a train test 3 kfold cross validation using the scikit-learn library in hopes of substantiating our derived R^2 value from the statsmodels.api exploration.

The following code cells below is the code used for these methods
"""

import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import itertools
!apt-get install libcairo2-dev libjpeg-dev libgif-dev
!pip install pycairo
import cairo
!pip install python-igraph
import igraph
import requests
!pip install scikit-learn
import statsmodels.api as sm
from sklearn.model_selection import KFold
from sklearn.preprocessing import PolynomialFeatures
from scipy.stats import t
import warnings

# google drive link should work for pulling this ~95mb data set. If link is broken
# please email me muellese@oregonstate.edu

url = 'https://drive.google.com/uc?id=1bj1ovhlwZ3ihwBD8BtSAxAZth8TEAcQy&export=download'

response = requests.get(url)

with open('Maize_PPI.txt', 'wb') as f:
  f.write(response.content)

ppi_df = pd.read_csv('Maize_PPI.txt', sep='\t', names=['protein1','protein2'])

print(ppi_df)

with open('Maize_mutation_info.txt', 'wb') as f:
  f.write(response.content)

mutation_df = pd.read_csv('Maize_mutation_info.txt', names=['allele','gene_id', 'fitness cost of mutation'])

print(mutation_df)

"""**Data wrangling**

Going to remove the protein assay tags "_P01" etc from the end of each entry. We are doing this so that the node names will match the names of our input gene dataset
"""

# clean tags from protein df
ppi_df['protein1'] = ppi_df['protein1'].str.split("_", n=1).str.get(0)
ppi_df['protein2'] = ppi_df['protein2'].str.split("_", n=1).str.get(0)

# there are some tags in mutation df that need to be cleaned too
mutation_df['gene_id'] = mutation_df['gene_id'].str.split("_", n=1).str.get(0)

"""**Make graph object and calculate betweeness**"""

g = igraph.Graph.TupleList(ppi_df.values.tolist(), directed=False)
print(g.summary())

"""**More Wrangling**

Make boolean mask to get filtered mutation df for genes only present in ppi
"""

gene_list = mutation_df['gene_id'].values.tolist()

vertex_list = []
for value in gene_list:
  if (ppi_df['protein1'] == value).any():
    vertex_list.append(value)
  elif (ppi_df['protein2'] == value).any():
    vertex_list.append(value)

g_betweeness = g.betweenness(vertices=vertex_list, directed=False) # calculate betweeness

"""**Print some betweeness centrality statistics**"""

print(f'The protein with the highest betweenness centrality is {g.vs["name"][np.argmax(g_betweeness)]} with a score of {np.max(g_betweeness):0.2f}')

"""**Lets calculate degrees and avg neighbor degrees**"""

degrees = [g.degree(node) for node in vertex_list]

avg_neighbor_degrees = []

for name in vertex_list:
    vertex = g.vs.find(name)  # Find the vertex object by name
    if vertex is not None:
        neighbors = g.neighbors(vertex.index)
        neighbor_degrees = [g.degree(n) for n in neighbors]
        avg_neighbor_degree = sum(neighbor_degrees) / len(neighbor_degrees) if neighbor_degrees else 0
        avg_neighbor_degrees.append((name, avg_neighbor_degree))

"""**Zip the names to the betweeness value**

We use this to add a betweeness column to the mutation df for statistical analysis.
"""

betweeness_dict = dict(zip(vertex_list,g_betweeness))

degree_dict = dict(zip(vertex_list, degrees))

avg_n_degrees_dict = dict(avg_neighbor_degrees)

"""**Build df and merge to mutation df to produce final working df for statisitical analysis**"""

genes = []
betweeness = []
for key, value in betweeness_dict.items():
  genes.append(key)
  betweeness.append(value)

betweeness_df = pd.DataFrame({'gene_id': genes, 'betweeness': betweeness})

final_df = pd.merge(mutation_df, betweeness_df, on = 'gene_id', how='inner')

genes = []
degrees = []
for key, value in degree_dict.items():
  genes.append(key)
  degrees.append(value)

degrees_df = pd.DataFrame({'gene_id': genes, 'degree': degrees})

final_df = pd.merge(final_df, degrees_df, on = 'gene_id', how='inner')

genes = []
avg_ndegrees = []
for key, value in avg_n_degrees_dict.items():
  genes.append(key)
  avg_ndegrees.append(value)

avg_ndegrees_df = pd.DataFrame({'gene_id': genes, 'avg neigbor degrees': avg_ndegrees})

final_df = pd.merge(final_df, avg_ndegrees_df, on = 'gene_id', how='inner')

"""## **Run statistical analysis on feature found: betweeness, avg neighbor degree, and degree**

**Multiple Regression**
"""

#try multiple regression model
X = sm.add_constant(final_df[['betweeness', 'degree', 'avg neigbor degrees']])
multiple_regression = sm.OLS(final_df['fitness cost of mutation'], X).fit()
print(multiple_regression.summary())

"""**Looks like only betweeness is significant**

Lets run polynomial regression to see if we can fit better
"""

for degree in [2,3,4]:
    # Create polynomial features
    X_poly = np.column_stack([final_df['betweeness']**i for i in range(1, degree + 1)])
    X_poly = sm.add_constant(X_poly)  # Add constant term for intercept

    # Fit polynomial regression
    model = sm.OLS(final_df['fitness cost of mutation'], X_poly).fit()

    # Prints
    print(f"\nPolynomial Regression (Degree {degree}):")
    print("R-squared:", model.rsquared)
    print("Adjusted R-squared:", model.rsquared_adj)
    print("AIC:", model.aic)
    print("BIC:", model.bic)

"""**Choosing 2nd degree polynomial**"""

X_poly = np.column_stack([final_df['betweeness']**i for i in range(1, 3)])
X_poly = sm.add_constant(X_poly)  # Add constant term for intercept

# Fit polynomial regression model
model = sm.OLS(final_df['fitness cost of mutation'], X_poly).fit()

print(model.summary())

"""**Plot second degree polynomial model**"""

y_pred = model.predict(X_poly)
x = final_df['betweeness']
y = final_df['fitness cost of mutation']
plt.scatter(x, y, label='data')

# x values for regression line
sort_indices = np.argsort(x)
plt.plot(x[sort_indices], y_pred[sort_indices], color='red', label='Second-degree Regression')

plt.xlabel('Betweeness')
plt.ylabel('Fitness Cost of Mutation')
plt.legend()
plt.title('Scatter Plot with Second-degree Polynomial Regression')
plt.show()

"""**Plot Observed vs. Predicted**"""

predictions = model.predict(X_poly)
y = final_df['fitness cost of mutation']
plt.scatter(final_df['fitness cost of mutation'], predictions)
plt.plot([y.min(), y.max()], [y.min(), y.max()], color='red', label='1:1 refrence line')
plt.xlabel('Observed')
plt.ylabel('Predicted')
plt.title('Observed vs Predicted')
plt.legend()
plt.show()

"""**Run Polynomial regression on train test split to explore prediction power of model**"""

num_folds = 3

# make kfold method and arrays for storage
kf = KFold(n_splits=num_folds)
train_r2_scores = []
test_r2_scores = []

for train_index, test_index in kf.split(X):
    # Split for train and test
    X_train, X_test = final_df['betweeness'].iloc[train_index], final_df['betweeness'].iloc[test_index]
    y_train, y_test = final_df['fitness cost of mutation'].iloc[train_index], final_df['fitness cost of mutation'].iloc[test_index]

    # build on fold
    X_train_poly = sm.add_constant(np.column_stack((X_train, X_train**2)))
    X_test_poly = sm.add_constant(np.column_stack((X_test, X_test**2)))

    # fit model
    train_r2 = sm.OLS(y_train, X_train_poly).fit().rsquared
    train_r2_scores.append(train_r2)

    # find r squared
    test_r2 = sm.OLS(y_test, X_test_poly).fit().rsquared_adj
    test_r2_scores.append(test_r2)

# find mean and std for fold arrays and print
mean_train_r2 = np.mean(train_r2_scores)
mean_test_r2 = np.mean(test_r2_scores)
std_train_r2 = np.std(train_r2_scores)
std_test_r2 = np.std(test_r2_scores)

print(f"Mean train R-squared: {mean_train_r2:.4f} +/- {std_train_r2:.4f}")
print(f"Mean test R-squared: {mean_test_r2:.4f} +/- {std_test_r2:.4f}")

"""**Run permutation test to explore distribution of p-values**

"""

warnings.filterwarnings("ignore")

def permutation_test(data, model, num_permutations):
    observed_coefficients = model.params
    num_features = len(observed_coefficients)

    # for significant returns
    num_extreme = np.zeros(num_features)

    # Fit second degree polynomial
    y = data['fitness cost of mutation']
    X_poly = np.column_stack([data['betweeness']**i for i in range(1, 3)])
    X_poly = sm.add_constant(X_poly)  # Add constant term for intercept

    for _ in range(num_permutations):
        np.random.shuffle(y)

        #fit for shuffle
        permuted_model = sm.OLS(y, X_poly).fit()
        permuted_coefficients = permuted_model.params

        # Check for more extreme than observed
        num_extreme += np.abs(permuted_coefficients) >= np.abs(observed_coefficients)

    p_values = (num_extreme + 1) / (num_permutations + 1)

    return p_values

data = final_df[['betweeness', 'fitness cost of mutation']]

y = data['fitness cost of mutation']
X_poly = np.column_stack([data['betweeness']**i for i in range(1, 3)])
X_poly = sm.add_constant(X_poly)
model = sm.OLS(y, X_poly).fit()

# Run permutation test
num_permutations = 10000
permutation_p_values = permutation_test(data, model, num_permutations)
results = pd.DataFrame({'permutation_p_values': permutation_p_values})

print(results)

"""## **Results**

After working with data, we ended up with three potential explanatory variables that might be associated with the changes in the fitness cost of mutation: betweenness centrality, degree, and average neighbor degree.
First, we run a multiple linear regression model with the fitness cost of mutation as a response variable and betweenness centrality, degree, and average neighbor degree as covariates to test the hypothesis that the effect of any of them is different from zero.

According to the multiple linear regression results, we have enough statistical evidence to reject the null hypothesis and assume that, on average, a one-unit increase in the betweenness centrality leads to an increase in the expected fitness cost of mutation by 2.105e-06 at the 5% significance level. At the same time, we do not have enough statistical evidence to assume an association between the fitness cost, on the one hand, and degree and average neighbor degrees, on the other hand.

Though the multiple linear regression yielded some statistically significant results, the current model explains only about 5% of the whole variation in the fitness cost. Given the limited number of significant covariates, one of the possible ways to deepen the analysis is by exploring polynomial models to see if the relationship between the betweenness centrality and the fitness cost is linear.

We created three additional models based on adding quadratic, cubic, and quartic terms. Adding more polynomials did not produce statistically significant results and increased the risk of overfitting. To choose the model, we looked at the change in the adjusted R-squared coefficient to see if a model explains more variation (1); AIC and BIC values to compare model fits (which is asymptotically equivalent to cross-validation) (2).

According to the results, the quadratic and cubic model looks the most optimal, explaining 11% and 19% of variation, respectively. However, the cubic model showed slightly higher AIC and BIC values, and if we are to approach the results conservatively, it is more reliable to stick to the quadratic model. The quadratic regression output showed statistically significant results at the 5% significance level, based on which we can assume that the possible impact of the betweenness centrality on the fitness cost is diminishing after a certain threshold.

Finally, we performed two validations. First, we  cross-validated the R-squared coefficients using the K-fold approach due to our data set's limited number of observations, which produced similar results in both train and test datasets. Secondly, we performed a permutation test for 10,000 iteration to verify that our signifcant coefficent p-values were not due to breakage of underlying assumptions for our parametric test. The permutaiton test was successful returning significant p values -- indicating that there was a significant amount of non-signifcant results from random shuffling. This means our p-values are still valid even if the underlying assumptions of the parametric test were broken. It is still essential to understand that despite being statistically significant, other covariates are contributing to the change in the fitness cost not captured by our model, and further investigation is required.

## **Conclusion**

This study found a significant association between betweeness centrality of the protein and cooresponing fitness cost of a mutation in the associated protein coding gene. While the R^2 value of the cooresponding regression model was low this does not negate the fact that this analysis elucidated a statisitcally significant effect between betweeness centrality of a protein and the fitness effect of a mutation of the associated gene.

These finding suggest that disruption of proteins that serve as "bridges" in the PPI of other proteins may be related to overal affect those genes have with maize fitness.

A drawback of this study was the lack of datapoints.  We did not have many  measurements of fitness and this may have limited the prediction capabilities of our model.

## Citations and Acknowledgments

All data and background information was sources from an unrealeased paper from the OSU Fowler Lab.

Thank Dr. John Fowler for generously allowing us to use his muation data.

PPI network citation:

G. Zhu, R. Jiang and X. -M. Zhao, "PPIM: A protein-protein interaction database for Maize," 2017 13th IEEE Conference on Automation Science and Engineering (CASE), Xi'an, China, 2017, pp. 97-97, doi: 10.1109/COASE.2017.8256085.
"""