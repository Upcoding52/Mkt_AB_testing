

[globals().pop(var) for var in dir() if not var.startswith("__")]

from math import sqrt,ceil,log
import pandas as pd
import scipy.stats


data = pd.read_csv('AB_test_data.csv')
data.drop(['date','id'], axis=1, inplace=True)
A = data[data.Variant == 'A']
A.drop(['Variant'], axis=1, inplace=True)
B = data[data.Variant == 'B']
B.drop(['Variant'], axis=1, inplace=True)

p_A = len(A[A.purchase_TF == True])/len(A)
p_B = len(B[B.purchase_TF == True])/len(B)
p = (p_A + p_B)/2

n_sample = 10

#2
alpha = 0.05
beta = 0.2
delta = p_B - p_A
t_alpha = scipy.stats.norm.ppf(alpha/2)
t_beta = scipy.stats.norm.ppf(beta)

n_star = ceil((t_alpha*sqrt(2*p*(1-p))+t_beta*sqrt(p_A*(1-p_A)+p_B*(1-p_B)))**2*(1/delta**2))
print('Optimal sample size is %f'%n_star)  

# Boundaries
up = log(1/alpha)
low = log(beta)

count_all = 0
for i in range(n_sample):    
    A_sample = A.sample(n = n_star, random_state=i)
    B_sample = B.sample(n = n_star, random_state=i)
    A_sample_p = len(A_sample[A_sample.purchase_TF == True])/len(A_sample)
    B_sample_p = len(B_sample[B_sample.purchase_TF == True])/len(B_sample)
    t = scipy.stats.norm.sf(abs(A_sample_p-B_sample_p)/sqrt((p_A*(1-p_A)+p_B*(1-p_B))/n_star))*2
    print('type1 test p-value is %f'%t)
    print('type2 test p-value is %f'%scipy.stats.norm.sf(abs(B_sample_p-A_sample_p-delta)
                                                         /sqrt((p_A*(1-p_A)+p_B*(1-p_B))/n_star)))
    lamb = 0
    count = 0
    
    for j in B_sample['purchase_TF']:
         
        if lamb > low and lamb < up:
            if j == True:
                lamb += log(p_B/p_A)
                count += 1
            else:
                lamb += log((1-p_B)/(1-p_A))
                count += 1
        if lamb <= low:
            print('In sample %f, we accept H0'%(i+1))
            break
        if lamb >= up:
            print('In sample %f, we accept H1'%(i+1))
            break
    print('The stop step in sample %f is %f'%(i+1,count))
    print('-------------------------------------------------------------')
    count_all += count

print('Average stop steps is %f'%(count_all/n_sample)) 
      
        
        
