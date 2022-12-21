import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('rmse_data.csv')

plt.figure()
for preconditioner in ['id', 'diag', 'ichol']:
    for coarse_init in [True, False]:
        df_filtered = df.loc[df['preconditioner']==preconditioner].loc[df['coarseInit']==coarse_init]
        max_iter = df_filtered['maxIter']
        sample_rmse = df_filtered['rmse_sample']
        plt.loglog(max_iter, sample_rmse, label=r'$P = {}, u_0 = {}$'.format(preconditioner, 'u_c' if coarse_init else 'None'))
        plt.xticks(max_iter, max_iter)
plt.title('RMSE of the first drawn sample')
plt.xlabel('number of iterations')
plt.ylabel('Single posterior sample RMSE')
plt.minorticks_off()
plt.legend()
plt.savefig('img/sample-rmse-plot.png', dpi=450)
plt.show()

plt.figure()
for preconditioner in ['id', 'diag', 'ichol']:
    for coarse_init in [True, False]:
        df_filtered = df.loc[df['preconditioner']==preconditioner].loc[df['coarseInit']==coarse_init]
        max_iter = df_filtered['maxIter']
        sample_rmse = df_filtered['rmse_std']
        plt.loglog(max_iter, sample_rmse, label=r'$P = {}, u_0 = {}$'.format(preconditioner, 'u_c' if coarse_init else 'None'))
        plt.xticks(max_iter, max_iter)
plt.title('RMSE of the posterior sample std')
plt.xlabel('number of iterations')
plt.ylabel('Posterior sample std RMSE')
plt.minorticks_off()
plt.legend()
plt.savefig('img/std-rmse-plot.png', dpi=450)
plt.show()
