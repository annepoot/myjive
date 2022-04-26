import numpy as np

from names import GPActions as gpact
from names import GPParamNames as gppn
from names import GlobNames as gn
from names import PropNames as prn
from model import Model, ModelFactory

NOBS = 'nobs'
OBSNOISE = 'obsNoise'
ALPHA = 'alpha'
RANDOMOBS = 'randomObs'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
DOFTYPES = ['dx']

class GPModel(Model):

    def take_action(self, action, params, globdat):
        print('GPModel taking action', action)

        if action == gpact.GETPRIORMEAN:
            self._get_prior_mean(params, globdat)
        elif action == gpact.GETPRIORCOVARIANCE:
            self._get_prior_covariance(params, globdat)
        elif action == gpact.GETPOSTERIORMEAN:
            self._get_posterior_mean(params, globdat)
        elif action == gpact.GETPOSTERIORCOVARIANCE:
            self._get_posterior_covariance(params, globdat)
        elif action == gpact.GETPRIORSAMPLES:
            self._get_prior_samples(params, globdat)
        elif action == gpact.GETPOSTERIORSAMPLES:
            self._get_posterior_samples(params, globdat)

    def configure(self, props, globdat):

        self._dc = globdat[gn.DOFSPACE].dof_count()

        # Get the number of observations either as a ratio or as a fixed integer
        self._nobs = float(props.get(NOBS,self._dc))
        if self._nobs < 1:
            self._nobs = self._dc * self._nobs
        self._nobs = int(np.round(self._nobs))

        # Get the observational noise
        self._noise2 = float(props.get(OBSNOISE))**2

        # Get the alpha parameter, or use the optimal alpha
        self._alpha = props.get(ALPHA, 'opt')
        if self._alpha != 'opt':
            self._alpha = float(self._alpha)

        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][prn.TYPE], props[SHAPE][INTSCHEME])
        self._rank = 1

    def configure_fem(self, globdat):

        # Get K, M and f from globdat
        K = globdat.get(gn.MATRIX0)
        M = globdat.get(gn.MATRIX2)
        f = globdat.get(gn.EXTFORCE)
        c = globdat.get(gn.CONSTRAINTS)

        # Constrain K, M and f
        self._Mc = c.constrain(M, f)[0]
        self._Kc, self._fc = c.constrain(K, f)

        # Get phi from Globdat
        self._phi = self._get_phis()[0]

        self._phi = self._get_phi_lumped(globdat)
        globdat['phi'] = self._phi

        # Compute the observation vector
        self._y = self._phi.T @ self._fc

        # Set alpha to the optimal value if necessary
        if self._alpha == 'opt':
            self._alpha = self._get_alpha_opt()
            print('Setting alpha to optimal value: {:.4f}'.format(self._alpha))


    def _get_prior_mean(self, params, globdat):

        #########################
        # PRIOR MEAN ON u and f #
        #########################

        # Only a 0-mean prior has been implemented
        params[gppn.PRIORMEAN] = np.zeros(self._dc)


    def _get_posterior_mean(self, params, globdat):

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._alpha**2 * self._phi.T @ self._Mc @ self._phi + np.identity(self._nobs) * self._noise2

        # Do the cholesky decomposition of Sigma_obs only if necessary
        if not '_sqrtObs' in vars(self):
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

        # Solve the system for the observed forces
        if not '_v0' in vars(self):
            self._v0 = np.linalg.solve(self._sqrtObs, self._y)

        # Get the posterior of the force field
        if not '_f_post' in vars(self):
            self._f_post = self._alpha**2 * self._Mc @ self._phi @ np.linalg.solve(self._sqrtObs.T, self._v0)

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')

        if field == 'u':

            #######################
            # POSTERIOR MEAN ON u #
            #######################

            # u_bar = alpha2 * inv(K) * M * phi * inv(alpha2 * phi.T * M * phi.T + Sigma_e) * f_obs

            if not '_u_post' in vars(self):

                # Solve to get the posterior of the displacement field
                self._u_post = np.linalg.solve(self._Kc, self._f_post)

            # Return the posterior of the displacement field
            params[gppn.POSTERIORMEAN] = self._u_post

        elif field == 'f':

            #######################
            # POSTERIOR MEAN ON f #
            #######################

            # f_bar = alpha2 * M * phi * inv(alpha2 * phi.T * M * phi.T + Sigma_e) * f_obs

            # Return the posterior of the force field
            params[gppn.POSTERIORMEAN] = self._f_post

        else:

            raise ValueError(field)


    def _get_prior_covariance(self, params, globdat):

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')
        fullSigma = params.get(gppn.FULLCOVARIANCE, False)

        if field == 'u':

            #########################
            # PRIOR COVARIANCE ON u #
            #########################

            # Sigma = alpha2 * inv(K) * M * inv(K)

            # Do the cholesky decomposition of M only if necessary
            if not '_sqrtM' in vars(self):
                self._sqrtM = np.linalg.cholesky(self._Mc)

            # Solve the system for each dof
            if not '_V2' in vars(self):
                self._V2 = np.linalg.solve(self._Kc, self._sqrtM)

            # Check if the full covariance matrix should be returned
            if fullSigma:

                # If so, compute the full covariance matrix
                Sigma_prior = self._alpha**2 * self._V2 @ self._V2.T

            else:

                # Compute only the diagonal of the covariance matrix
                sigma_prior = np.zeros(self._dc)
                for i in range(self._dc):
                    sigma_prior[i] = self._alpha**2 * self._V2[i,:] @ self._V2[i,:].T

        elif field == 'f':

            #########################
            # PRIOR COVARIANCE ON f #
            #########################

            # Sigma = alpha2 * M

            # Check if the full covariance matrix should be returned
            if fullSigma:

                # If so, compute the full covariance matrix
                Sigma_prior = self._alpha**2 * self._Mc

            else:

                # If not, compute only the diagonal of the covariance matrix
                sigma_prior = self._alpha**2 * self._Mc.diagonal()

        else:

            raise ValueError(field)

        # Return either the full covariance matrix, or only its diagonal
        if fullSigma:
            params[gppn.PRIORCOVARIANCE] = Sigma_prior
        else:
            params[gppn.PRIORCOVARIANCE] = sigma_prior


    def _get_posterior_covariance(self, params, globdat):

        # Check if the posterior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')
        fullSigma = params.get(gppn.FULLCOVARIANCE, False)

        # Check if the prior variance is given in the params, otherwise, take an action to obtain it
        if not gppn.PRIORCOVARIANCE in params:
            self.take_action(gpact.GETPRIORCOVARIANCE, params, globdat)

        # Get the prior variance from the params
        sigma_prior = params[gppn.PRIORCOVARIANCE]

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._alpha**2 * self._phi.T @ self._Mc @ self._phi + np.identity(self._nobs) * self._noise2

        # Do the cholesky decomposition of Sigma_obs only if necessary
        if not '_sqrtObs' in vars(self):
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

        # Solve the inverse observation covariance once for each observation
        if not '_V1' in vars(self):
            self._V1 = np.linalg.solve(self._sqrtObs, self._phi.T @ self._Mc)

        if field == 'u':

            #############################
            # POSTERIOR COVARIANCE ON u #
            #############################

            # Sigma = alpha2 * inv(K) * M * inv(K) - alpha4 * inv(K) * M * phi * inv(alpha2 * phi.T * M * phi.T + Sigma_e) * phi.T * M * inv(K)

            # Solve the system for each dof
            if not '_V3' in vars(self):
                self._V3 = np.linalg.solve(self._Kc, self._V1.T)

            # Check if the full covariance matrix should be returned
            if fullSigma:

                # If so, compute the full covariance matrix
                Sigma_post = sigma_prior.copy()
                Sigma_post -= self._alpha**4 * self._V3 @ self._V3.T

            else:

                # If not, compute only the diagonal of the covariance matrix
                sigma_post = sigma_prior.copy()
                for i in range(self._dc):
                    sigma_post[i] -= self._alpha**4 * self._V3[i,:] @ self._V3[i,:].T

        elif field == 'f':

            #############################
            # POSTERIOR COVARIANCE ON f #
            #############################

            # Sigma = alpha2 * M - alpha4 * M * phi * inv(alpha2 * phi.T * M * phi + Sigma_e) * phi.T * M

            # Check if the full covariance matrix should be returned
            if fullSigma:

                # If so, compute the full covariance matrix
                Sigma_post = sigma_prior.copy()
                Sigma_post -= self._alpha**4 * self._V1.T @ self._V1

            else:

                # If not, compute only the diagonal of the covariance matrix
                sigma_post = sigma_prior.copy()
                for i in range(self._dc):
                    sigma_post[i] -= self._alpha**4 * self._V1[:,i].T @ self._V1[:,i]

        else:

            raise ValueError(field)

        # Return either the full covariance matrix, or only its diagonal
        if fullSigma:
            params[gppn.POSTERIORCOVARIANCE] = Sigma_post
        else:
            params[gppn.POSTERIORCOVARIANCE] = sigma_post


    def _get_prior_samples(self, params, globdat):

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')
        nsamples = params.get(gppn.NSAMPLE, 1)

        # Do the cholesky decomposition of M only if necessary
        if not '_sqrtM' in vars(self):
            self._sqrtM = np.linalg.cholesky(self._Mc)

        # Define the array that stores the samples
        samples = np.zeros((self._dc, nsamples))

        # Get a loop to create all samples
        for i in range(nsamples):

            # Get a sample from the standard normal distribution
            z = np.random.randn(self._dc)

            # Get the sample of the force field
            f = self._alpha * self._sqrtM @ z

            if field == 'u':

                #####################
                # PRIOR SAMPLE ON u #
                #####################

                # Compute the corresponding displacement field
                u = np.linalg.solve(self._Kc, f)

                samples[:,i] = u

            elif field == 'f':

                #####################
                # PRIOR SAMPLE ON f #
                #####################

                samples[:,i] = f

            else:

                raise ValueError(field)

        # Return the array of samples
        params[gppn.PRIORSAMPLES] = samples


    def _get_posterior_samples(self, params, globdat):

        # Check if the posterior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')
        nsamples = params.get(gppn.NSAMPLE, 1)

        # Do the cholesky decomposition of M only if necessary
        if not '_sqrtM' in vars(self):
            self._sqrtM = np.linalg.cholesky(self._Mc)

        # Do the cholesky decomposition of M only if necessary
        if not '_sqrtM' in vars(self):
            self._sqrtM = np.linalg.cholesky(self._Mc)

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._alpha**2 * self._phi.T @ self._Mc @ self._phi + np.identity(self._nobs) * self._noise2

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
            z1 = np.random.randn(self._dc)
            z2 = np.random.randn(self._nobs)

            # Get x1 ~ N(0, alpha^2 M) and x2 ~ N(0, Sigma_e)
            x1 = self._alpha * self._sqrtM @ z1
            x2 = self._sqrtNoise @ z2

            # Compute the perturbed observation
            f_pert = self._y - self._phi.T @ x1 + x2

            # Multiply the perturbed observation by the Kalman gain
            f = np.linalg.solve(self._sqrtObs, f_pert)
            f = np.linalg.solve(self._sqrtObs.T, f)
            f = self._alpha**2 * self._Mc @ self._phi @ f

            # Add the prior sample
            f += x1

            if field == 'u':

                #########################
                # POSTERIOR SAMPLE ON u #
                #########################

                # Compute the corresponding displacement field
                u = np.linalg.solve(self._Kc, f)

                samples[:,i] = u

            elif field == 'f':

                #########################
                # POSTERIOR SAMPLE ON f #
                #########################

                samples[:,i] = f

            else:

                raise ValueError(field)

        # Return the array of samples
        params[gppn.POSTERIORSAMPLES] = samples


    def _get_alpha_opt(self):

        #################
        # OPTIMAL ALPHA #
        #################

        # If so, determine the optimal value of alpha
        L = np.linalg.cholesky(self._phi.T @ self._Mc @ self._phi)
        v = np.linalg.solve(L, self._y)
        alpha2 = v.T @ v / self._nobs

        return np.sqrt(alpha2)


    def _get_phis(self):

        # Store the number of dofs, and number of observatiosn
        dc = self._dc
        nobs = self._nobs

        # Create an empty array
        phi = np.zeros((dc, nobs))
        phi_sub = np.zeros((dc, nobs))

        # Get the indices of the nodes that will be observed
        vec_obs = np.rint(np.linspace(0, dc-1, nobs))
        vec_obs = vec_obs.astype(int)

        for i in range(nobs):

            phi[vec_obs[i], i] = 1
            phi_sub[vec_obs[i], i] = 1

        j = 0
        j_old = 0

        for i in range(nobs-1):
            while j < dc:
                if np.isclose(phi[j,i+1], 1):
                    phi[j_old:j+1,i] = np.linspace(1, 0, j-j_old+1)
                    phi[j_old:j+1,i+1] = np.linspace(0, 1, j-j_old+1)
                    j_old = j
                    break

                j = j+1

        return phi, phi_sub


    def _get_phi_lumped(self, globdat):

        suffix = 'Coarse'

        elemsc = globdat[gn.ESET + suffix]
        elems = globdat[gn.ESET]
        nodesc = globdat[gn.NSET + suffix]
        nodes = globdat[gn.NSET]
        dofsc = globdat[gn.DOFSPACE + suffix]
        dofs = globdat[gn.DOFSPACE]

        phi = np.zeros((len(nodes), len(nodesc)))

        # Go over the coarse mesh
        for elemc in elemsc:
            inodesc = elemc.get_nodes()
            idofsc = dofsc.get_dofs(inodesc, DOFTYPES[0:self._rank])
            coordsc = np.stack([nodesc[i].get_coords() for i in inodesc], axis=1)[0:self._rank, :]

            lower = np.min(coordsc)
            upper = np.max(coordsc)

            # Go over the fine mesh
            for elem in elems:
                inodes = elem.get_nodes()
                idofs = dofs.get_dofs(inodes, DOFTYPES[0:self._rank])
                coords = np.stack([nodes[i].get_coords() for i in inodes], axis=1)[0:self._rank, :]

                # Check if the elements overlap
                if np.max(coords) > lower and np.min(coords) < upper:

                    # Go over the nodes of the fine element
                    for n in range(len(inodes)):

                        # Get the nodal coords
                        coord = coords[:, n][0]

                        # Check if the node overlaps
                        if coord >= lower and coord <= upper:

                            # If so, get the relative position of the node
                            relcoord = -1 + 2 * (coord - lower) / (upper - lower)

                            # Get the shape function values at the location of the coords
                            svals = self._shape.eval_shape_values(relcoord)

                            # Get the dofs belonging to this node
                            idof = dofs.get_dofs([inodes[n]], DOFTYPES[0:self._rank])

                            # Store the relative shape function values in the phi matrix
                            phi[idof, idofsc] = svals

        return phi

def declare(factory):
    factory.declare_model('GP', GPModel)
