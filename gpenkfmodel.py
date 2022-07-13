import numpy as np

from names import Actions as act
from names import GPActions as gpact
from names import ParamNames as pn
from names import GPParamNames as gppn
from names import GlobNames as gn
from names import PropNames as prn
from gpmodel import GPModel

NOBS = 'nobs'
OBSNOISE = 'obsNoise'
ALPHA = 'alpha'
BETA = 'beta'
PRIOR = 'prior'
TYPE = 'type'
FUNC = 'func'
HYPERPARAMS = 'hyperparams'
RANDOMOBS = 'randomObs'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
PDNOISE = 'pdNoise'
ENSEMBLE = 'ensemble'


class GPEnKFModel(GPModel):

    def configure(self, props, globdat):

        super().configure(props, globdat)

        self._nens = props.get(ENSEMBLE, 100)


    def _configure_prior(self, params, globdat):

        # !!! Note: here, Sigma is still being assembled.
        # How to get the ensemble without sampling from the prior?
        super()._configure_prior(params, globdat)

        # Get the ensemble from the prior distribution
        params = {gppn.NSAMPLE:self._nens}

        super()._get_prior_samples(params, globdat)

        # Get the ensemble back from the array
        self._X = params[gppn.PRIORSAMPLES]

        # Get the mean of X, and the deviation from the mean
        EX = np.mean(self._X, axis=1)
        self._A = self._X - np.tile(EX, (self._nens,1)).T

        # Compute the product of H and A (which is nobs x nens)
        self._HA = self._H @ self._A

        # Delete Sigma altogether, to make sure it's never used
        del self._Sigma


    def _get_prior_mean(self, params, globdat):

        # Return the prior of the displacement field
        params[gppn.PRIORMEAN] = self._m


    def _get_posterior_mean(self, params, globdat):

        ##################
        # POSTERIOR MEAN #
        ##################

        # u_bar = Sigma * H.T * inv(H * Sigma * H.T + Sigma_e) * y

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._HA @ self._HA.T / (self._nens-1) + np.identity(self._nobs) * self._noise2

        # Do the cholesky decomposition of Sigma_obs only if necessary
        if not '_sqrtObs' in vars(self):
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

        # Solve the system for the observed forces
        if not '_v0' in vars(self):
            self._v0 = np.linalg.solve(self._sqrtObs, self._y)

        # Get the posterior of the force field
        if not '_u_post' in vars(self):
            self._u_post = self._m + self._A @ self._HA.T @ np.linalg.solve(self._sqrtObs.T, self._v0) / (self._nens-1)

        # Return the posterior of the displacement field
        params[gppn.POSTERIORMEAN] = self._u_post


    def _get_prior_covariance(self, params, globdat):

        fullSigma = params.get(gppn.FULLCOVARIANCE, False)

        ####################
        # PRIOR COVARIANCE #
        ####################

        # Check if the full covariance matrix should be returned
        if fullSigma:

            # If so, compute the full covariance matrix
            Sigma_prior = self._A @ self._A.T / (self._nens - 1)
            params[gppn.PRIORCOVARIANCE] = Sigma_prior

        else:

            # If not, compute only the diagonal of the covariance matrix
            var_prior = np.zeros(self._dc)
            for i in range(self._dc):
                var_prior[i] = self._A[i,:] @ self._A[i,:].T / (self._nens - 1)

            params[gppn.PRIORCOVARIANCE] = var_prior


    def _get_posterior_covariance(self, params, globdat):

        fullSigma = params.get(gppn.FULLCOVARIANCE, False)

        ########################
        # POSTERIOR COVARIANCE #
        ########################

        # Sigma_bar = Sigma - Sigma * H.T * inv(H * Sigma * H.T + Sigma_e) * H * Sigma

        # Check if the prior variance is given in the params, otherwise, take an action to obtain it
        if not gppn.PRIORCOVARIANCE in params:
            self.take_action(gpact.GETPRIORCOVARIANCE, params, globdat)

        # Get the prior variance from the params
        var_prior = params[gppn.PRIORCOVARIANCE]

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._HA @ self._HA.T / (self._nens-1) + np.identity(self._nobs) * self._noise2

        # Do the cholesky decomposition of Sigma_obs only if necessary
        if not '_sqrtObs' in vars(self):
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

        # Solve the inverse observation covariance once for each observation
        if not '_V1' in vars(self):
            self._V1 = np.linalg.solve(self._sqrtObs, self._HA @ self._A.T / (self._nens-1))

        # Check if the full covariance matrix should be returned
        if fullSigma:

            # If so, compute the full covariance matrix
            Sigma_post = var_prior.copy()
            Sigma_post -= self._V1.T @ self._V1
            params[gppn.POSTERIORCOVARIANCE] = Sigma_post

        else:

            # If not, compute only the diagonal of the covariance matrix
            var_post = var_prior.copy()
            for i in range(self._dc):
                var_post[i] -= self._V1[:,i].T @ self._V1[:,i]

            params[gppn.POSTERIORCOVARIANCE] = var_post


    def _get_prior_samples(self, params, globdat):

        nsamples = params.get(gppn.NSAMPLE, 1)
        rng = params.get(gppn.RNG, np.random.default_rng())

        #################
        # PRIOR SAMPLES #
        #################

        # u = m + sqrt(Sigma) * z

        # Do the cholesky decomposition of Sigma only if necessary
        if not '_sqrtSig' in vars(self):
            self._sqrtSig = self._A / np.sqrt(self._nens-1)

        # Define the array that stores the samples
        samples = np.zeros((self._dc, nsamples))

        # Get a loop to create all samples
        for i in range(nsamples):

            # Get a sample from the standard normal distribution
            z = rng.standard_normal(self._dc)

            # Get the sample of the force field
            u = self._sqrtSig @ z + self._m

            samples[:,i] = u

        # Return the array of samples
        params[gppn.PRIORSAMPLES] = samples


    def _get_posterior_samples(self, params, globdat):

        nsamples = params.get(gppn.NSAMPLE, 1)
        rng = params.get(gppn.RNG, np.random.default_rng())

        #####################
        # POSTERIOR SAMPLES #
        #####################

        # u = m + sqrt(Sigma) * z1 + Sigma * H.T * inv(H * Sigma * H + Sigma_e) * (y - H * sqrt(Sigma) * z1 + sqrt(Sigma_e) * z2)

        # Do the cholesky decomposition of Sigma only if necessary
        if not '_sqrtSig' in vars(self):
            self._sqrtSig = self._A / np.sqrt(self._nens-1)

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._HA @ self._HA.T / (self._nens-1) + np.identity(self._nobs) * self._noise2

        # Do the cholesky decomposition of Sigma_obs only if necessary
        if not '_sqrtObs' in vars(self):
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

        # Do the cholesky decomposion of the noise only if necesary
        if not '_sqrtNoise' in vars(self):
            self._sqrtNoise = np.sqrt(self._noise2) * np.identity(self._nobs)

        # Define the array that stores the samples
        samples = np.zeros((self._dc, nsamples))

        # Get a loop to create all samples
        for i in range(nsamples):

            # Get two samples from the standard normal distribution
            z1 = rng.standard_normal(self._dc)
            z2 = rng.standard_normal(self._nobs)

            # Get x1 ~ N(0, alpha^2 M) and x2 ~ N(0, Sigma_e)
            x1 = self._sqrtSig @ z1
            x2 = self._sqrtNoise @ z2

            # Compute the perturbed observation
            u_pert = self._y - self._H @ x1 + x2

            # Multiply the perturbed observation by the Kalman gain
            u = np.linalg.solve(self._sqrtObs, u_pert)
            u = np.linalg.solve(self._sqrtObs.T, u)
            u = self._A @ self._HA.T @ u / (self._nens-1)

            # Add the prior sample and mean
            u += x1 + self._m

            samples[:,i] = u

        # Return the array of samples
        params[gppn.POSTERIORSAMPLES] = samples


    def _get_log_likelihood(self, params, globdat):

        ##################
        # LOG LIKELIHOOD #
        ##################

        # l = - 1/2 * y * inv(H * Sigma * H.T + Sigma_e) * y - 1/2 * log|H * Sigma * H.T + Sigma_e| - n/2 * log(2*pi)

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._HA @ self._HA.T / (self._nens-1) + np.identity(self._nobs) * self._noise2

        # Do the cholesky decomposition of Sigma_obs only if necessary
        if not '_sqrtObs' in vars(self):
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

        # Solve the system for the observed forces
        if not '_v0' in vars(self):
            self._v0 = np.linalg.solve(self._sqrtObs, self._y)

        l = - 0.5 * self._v0.T @ self._v0 - np.sum(np.log(self._sqrtObs.diagonal())) - 0.5 * self._nobs * np.log(2*np.pi)

        params[gppn.LOGLIKELIHOOD] = l


    def _get_phi(self, globdat):

        elems = globdat[gn.ESET]
        elemsc = globdat[gn.COARSEMESH][gn.ESET]
        nodes = globdat[gn.NSET]
        nodesc = globdat[gn.COARSEMESH][gn.NSET]
        dofs = globdat[gn.DOFSPACE]
        dofsc = globdat[gn.COARSEMESH][gn.DOFSPACE]

        phi = np.zeros((dofs.dof_count(), dofsc.dof_count()))

        # Go over the coarse mesh
        for elemc in elemsc:
            inodesc = elemc.get_nodes()
            idofsc = dofsc.get_dofs(inodesc, self._dof_types)
            coordsc = np.stack([nodesc[i].get_coords() for i in inodesc], axis=1)[0:self._rank, :]

            # Get the bounding box of the coarse element
            bbox = np.zeros((self._rank, 2))
            for i in range(self._rank):
                bbox[i,0] = min(coordsc[i,:])
                bbox[i,1] = max(coordsc[i,:])

            # Go over the fine mesh
            for elem in elems:
                inodes = elem.get_nodes()
                coords = np.stack([nodes[i].get_coords() for i in inodes], axis=1)[0:self._rank, :]

                # Check if the bounding boxes of the coarse and fine element overlap
                inside = True
                for i in range(self._rank):
                    if max(coords[i,:]) < bbox[i,0] or min(coords[i,:]) > bbox[i,1]:
                        inside = False
                        break

                # Check if the elements overlap
                if inside:

                    # Go over the nodes of the fine element
                    for n in range(len(inodes)):

                        # Get the nodal coords
                        coord = coords[:, n]

                        # Check if the node falls inside the bounding box
                        inside = True
                        for i in range(self._rank):
                            if coord[i] < bbox[i,0] or coord[i] > bbox[i,1]:
                                inside = False
                                break

                        if inside:

                            # Get the relative position of the node
                            loc_point = self._shape.get_local_point(coord, coordsc)

                            # Check if the node actually falls within shape
                            inside = self._shape.contains_local_point(loc_point, tol=1e-8)

                        # Only continue if both checks are passed
                        if inside:

                            # Get the shape function values at the location of the coords
                            svals = self._shape.eval_shape_functions(loc_point)

                            # Get the dofs belonging to this node
                            idofs = dofs.get_dofs([inodes[n]], self._dof_types)

                            # Store the relative shape function values in the phi matrix
                            for i, idof in enumerate(idofs):
                                phi[idof, idofsc[i::len(self._dof_types)]] = svals

        return phi


    def _get_variances(self, params, globdat):

        table = params[pn.TABLE]
        name = params[pn.TABLENAME]
        tbwts = params[pn.TABLEWEIGHTS]

        # Initialize the tables
        for dof_type in self._dof_types:
            if dof_type not in table:
                table[dof_type] = np.zeros(len(globdat[gn.NSET]))

        # Define a dictionary for the settings of u
        u_params = {gppn.FULLCOVARIANCE:False}

        if name == 'var_u_prior':
            self._get_prior_covariance(u_params, globdat)
            var = u_params[gppn.PRIORCOVARIANCE]
        elif name == 'var_u_post':
            self._get_posterior_covariance(u_params, globdat)
            var = u_params[gppn.POSTERIORCOVARIANCE]

        c = globdat[gn.CONSTRAINTS]
        cdofs, cvals = c.get_constraints()

        for elem in globdat[gn.ESET]:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, self._dof_types)

            for inode in inodes:
                idofs = globdat[gn.DOFSPACE].get_dofs([inode], self._dof_types)

                for i, dof_type in enumerate(self._dof_types):
                    table[dof_type][inode] = var[idofs[i]]

                tbwts[inode] = 1


    def _get_standard_deviations(self, params, globdat):
        table = params[pn.TABLE]
        name = params[pn.TABLENAME]

        params[pn.TABLENAME] = name.replace('std_', 'var_')

        # Get the variances first
        self._get_variances(params, globdat)

        for comp in table:
            table[comp] = np.sqrt(table[comp])

        params[pn.TABLENAME] = name


def declare(factory):
    factory.declare_model('GPEnKF', GPEnKFModel)
