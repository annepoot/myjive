import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla

from jive.solver.jit.cholesky import sparse_cholesky
from jive.solver.jit.spsolve import solve_triangular
from jive.fem.names import GPActions as gpact
from jive.fem.names import GPParamNames as gppn
from jive.fem.names import GlobNames as gn
import jive.util.proputils as pu

from gpmodel import GPModel

PRIOR = 'prior'
BOUNDARY = 'boundary'
GROUPS = 'groups'
DOFS = 'dofs'
COVS = 'covs'
DIAGONALIZED = 'diagonalized'
TYPE = 'type'

class GPfModel(GPModel):

    def configure(self, props, globdat):

        super().configure(props, globdat)

        # Get the type of BC enforcement
        bcprops = props.get(BOUNDARY, {})
        self._bctype = bcprops.get(TYPE, 'dirichlet')
        if self._bctype not in ['direct', 'dirichlet']:
            raise ValueError('boundary has to be "dirichlet" or "direct"')
        self._bcgroups = pu.parse_list(bcprops.get(GROUPS,'[]'))
        self._bcdofs = pu.parse_list(bcprops.get(DOFS, '[]'))
        self._bccovs = pu.parse_list(bcprops.get(COVS,'[]'), float)

        self._diagonalized = props[PRIOR].get(DIAGONALIZED, False)

    def _configure_fem(self, params, globdat):

        super()._configure_fem(params, globdat)

    def _configure_prior(self, params, globdat):

        # Get the observation operator
        self._Phic = spsp.csr_array(self._Phic)
        self._Phi = spsp.csr_array(self._Phi)
        self._H = self._Phic.T
        self._H = self._H.tocsr()

        # Check if alpha or beta should be optimized
        if len(self._hyperparams) == 1:
            key, value = list(self._hyperparams.items())[0]

            if value == 'opt':

                # Check if the only hyperparameter to optimize is a scaling in the front
                if self._covariance.startswith(key+'**2'):

                    # Get the covariance matrix without the hyperparameter
                    eval_dict = self._get_eval_dict(globdat)
                    Sigma = eval(self._covariance.replace(key+'**2', '1'), eval_dict)

                else:
                    raise ValueError('cannot find optimal value for ' + key)

                # Apply boundary conditions to the prior
                m, Sigma = self._apply_covariance_bcs(Sigma, globdat)

                self._hyperparams[key] = self._get_param_opt(m, Sigma)

        super()._configure_prior(params, globdat)

    def _get_prior_mean(self, params, globdat):

        # Get the prior mean on f
        super()._get_prior_mean(params, globdat)
        params[gppn.FIELD] = gn.EXTFORCE

    def _get_posterior_mean(self, params, globdat):

        # Get the posterior mean on f
        super()._get_posterior_mean(params, globdat)
        params[gppn.FIELD] = gn.EXTFORCE

    def _get_prior_covariance(self, params, globdat):

        # Get the prior covariance on f
        super()._get_prior_covariance(params, globdat)
        params[gppn.FIELD] = gn.EXTFORCE

    def _get_posterior_covariance(self, params, globdat):

        # Get the posterior covariance on f
        super()._get_posterior_covariance(params, globdat)
        params[gppn.FIELD] = gn.EXTFORCE

    def _get_prior_samples(self, params, globdat):

        # Get the prior samples on f
        super()._get_prior_samples(params, globdat)
        params[gppn.FIELD] = gn.EXTFORCE

    def _get_posterior_samples(self, params, globdat):

        # Get the posterior samples on f
        super()._get_posterior_samples(params, globdat)
        params[gppn.FIELD] = gn.EXTFORCE

    def _kalman_update(self, params, globdat):

        # Perform the Kalman update on f
        super()._kalman_update(params, globdat)
        params[gppn.FIELD] = gn.EXTFORCE

    def _get_param_opt(self, m_fc, Sigma_fc):

        # Determine the optimal value of alpha
        y = self._g - self._Phic.T @ m_fc
        L = sparse_cholesky(self._H @ Sigma_fc @ self._H.T + spsp.identity(self._nobs) * self._noise2)
        v = self._solve_triangular(L, y, lower=True)
        alpha2 = v.T @ v / self._nobs

        return np.sqrt(alpha2)

    def _get_eval_dict(self, globdat):

        # Define a dictionary with relevant functions
        eval_dict = {'inv':spspla.inv, 'exp':np.exp, 'norm':np.linalg.norm, 'np':np}
        eval_dict.update(self._hyperparams)

        # Check if we have an SPDE covariance
        if self._prior == 'SPDE':

            nodes = globdat[gn.NSET]
            dofspace = globdat[gn.DOFSPACE]

            if self._rank >= 1:
                x = np.zeros(self._dc)
            if self._rank >= 2:
                y = np.zeros(self._dc)

            for i in range(len(nodes)):
                icoords = nodes[i].get_coords()
                idofs = dofspace.get_dofs([i], dofspace.get_types())

                if self._rank >= 1:
                    x[idofs] = icoords[0]
                if self._rank >= 2:
                    y[idofs] = icoords[1]

            if self._rank >= 1:
                eval_dict['x'] = x
            if self._rank >= 2:
                eval_dict['y'] = y

            g = self._Phi @ np.linalg.solve((self._Phi.T @ self._Phi).toarray(), self._g)

            # Add the mass and stiffness matrices to the dictionary
            eval_dict['M'] = self._Mc
            eval_dict['K'] = self._Kc
            eval_dict['F'] = spsp.csr_array(np.outer(self._f, self._f))
            eval_dict['G'] = spsp.csr_array(np.outer(g, g))
            eval_dict['Phi'] = self._Phi

        return eval_dict

    def _apply_covariance_bcs(self, Sigma, globdat):
        Sigmac = Sigma.copy()
        mc = np.zeros(self._dc)

        # Add a tiny noise to ensure Sigma is positive definite rather than semidefinite
        Sigmac += self._pdnoise2 * spsp.identity(self._dc)

        # Check if the boundary condition should be applied directly or via dirichlet BCs
        if self._bctype == 'dirichlet':

            # Split K along boundary and internal nodes
            K_ib = -self._K[:,self._cdofs]
            K_ib[self._cdofs] = spsp.identity(len(self._cdofs))

            # Update the prior mean by observing the displacement at the bcs
            if self._mean == 'dirichlet':
                mc += K_ib @ self._cvals

            # Decouple the bc covariance from the internal nodes
            Sigmac[self._cdofs,:] *= 0.0
            Sigmac[:,self._cdofs] *= 0.0

            # Get a matrix that defines the constraint equations_set_arrayXarray
            conmat = spsp.lil_array((self._dc, len(self._bcgroups)))
            ds = globdat[gn.DOFSPACE]
            for i, (group, dof) in enumerate(zip(self._bcgroups, self._bcdofs)):
                idofs = ds.get_dofs(globdat[gn.NGROUPS][group], [dof])
                conmat[idofs, i] = 1
            conmat = conmat[self._cdofs,:]

            # Get the boundary covariance matrix
            covmat = spsp.diags(self._bccovs)**2
            Sigma_bc = conmat @ covmat @ conmat.T

            # Recouple the internal nodes based on the boundary covariance matrix
            Sigmac += K_ib @ Sigma_bc @ K_ib.T

        elif self._bctype == 'direct':
            raise ValueError('With GPfModel, BCs cannot be applied directly')

        else:
            raise ValueError('boundary has to be "dirichlet" or "direct"')

        # Add separate boundary noise to ensure positive definiteness
        noisediag = np.zeros(self._dc)
        noisediag[self._cdofs] = self._bcnoise2
        Sigmac += spsp.diags(noisediag)

        return mc, Sigmac

    def _solve_triangular(self, A, b, lower):
        return solve_triangular(A.tocsr(), b, lower=lower)

    def _get_Sigma_obs(self):

        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._H @ self._Sigma @ self._H.T + spsp.identity(self._nobs) * self._noise2

    def _get_sqrtObs(self):

        if not '_sqrtObs' in vars(self):
            self._get_Sigma_obs()
            self._sqrtObs = sparse_cholesky(self._Sigma_obs)

    def _get_sqrtSigma(self):

        if not '_sqrtSigma' in vars(self):
            if self._diagonalized:
                self._sqrtSigma = spsp.diags(np.sqrt(self._Sigma.sum(axis=1)), format='csr')
            else:
                self._sqrtSigma = sparse_cholesky(self._Sigma)

    def _get_sqrtNoise(self):

        if not '_sqrtNoise' in vars(self):
            self._sqrtNoise = np.sqrt(self._noise2) * spsp.identity(self._nobs)

    def _get_P(self):

        if not '_P' in vars(self):
            self._P = self._Phi @ spspla.inv(self._Phi.T @ self._Phi) @ self._Phi.T


def declare(factory):
    factory.declare_model('GPf', GPfModel)
