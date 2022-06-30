import numpy as np

from names import Actions as act
from names import GPActions as gpact
from names import ParamNames as pn
from names import GPParamNames as gppn
from names import GlobNames as gn
from names import PropNames as prn
from model import Model

NOBS = 'nobs'
OBSNOISE = 'obsNoise'
ALPHA = 'alpha'
RANDOMOBS = 'randomObs'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
# DOFTYPES = ['dx', 'dy']

class GPModel(Model):

    def take_action(self, action, params, globdat):
        print('GPModel taking action', action)

        if action == gpact.CONFIGUREFEM:
            self._configure_fem(params, globdat)
        elif action == gpact.GETPRIORMEAN:
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
        elif action == gpact.GETLOGLIKELIHOOD:
            self._get_log_likelihood(params, globdat)
        elif action == act.GETTABLE:
            if 'var' in params[pn.TABLENAME]:
                self._get_variances(params, globdat)
            elif 'std' in params[pn.TABLENAME]:
                self._get_standard_deviations(params, globdat)

    def configure(self, props, globdat):
        self._dc = globdat[gn.DOFSPACE].dof_count()
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][prn.TYPE], props[SHAPE][INTSCHEME])
        self._rank = self._shape.global_rank()

        # Get the observational noise
        self._noise2 = float(props.get(OBSNOISE))**2

        # Get the alpha parameter, or use the optimal alpha
        self._alpha = props.get(ALPHA, 'opt')
        if self._alpha != 'opt':
            self._alpha = float(self._alpha)

        # Get the dofs of the fine mesh
        self._dof_types = globdat[gn.DOFSPACE].get_types()
        # Add the dofs to the coarse mesh
        nodes = np.unique([node for elem in globdat[gn.COARSEMESH][gn.ESET] for node in elem.get_nodes()])
        for doftype in self._dof_types:
            globdat[gn.COARSEMESH][gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.COARSEMESH][gn.DOFSPACE].add_dof(node, doftype)

        # Get the number of observations (which is the number of dofs in the coarse mesh)
        self._nobs = globdat[gn.COARSEMESH][gn.DOFSPACE].dof_count()

    def _configure_fem(self, params, globdat):

        # Get K, M and f from globdat
        K = globdat.get(gn.MATRIX0)
        M = globdat.get(gn.MATRIX2)
        f = globdat.get(gn.EXTFORCE)
        c = globdat.get(gn.CONSTRAINTS)

        cdofs, cvals = c.get_constraints()

        # The posterior mean force vector has to contain the Dirichlet BCs
        mf = np.zeros_like(f)
        Kc, mf = c.constrain(K, mf)

        # Get the actual constrained stiffness matrix and force vector
        Kc, fc = c.constrain(K, f)

        # Get the mass matrix, and remove the Dirichlet BCs
        Mc = M.copy()
        Mc[cdofs,:] = Mc[:,cdofs] = 0.0
        Mc += 1e-10 * np.identity(Mc.shape[0])

        # Store all constrained matrices and vectors
        self._Mc = Mc
        self._Kc = Kc
        self._fc = fc
        self._m = np.linalg.solve(Kc, mf)

        # Get the phi matrix, and constrain the Dirichlet BCs
        phi = self._get_phi(globdat)

        self._Phi = phi.copy()
        globdat['Phi'] = self._Phi

        for i in range(phi.shape[1]):
            for cdof in cdofs:
                if np.isclose(phi[cdof,i], 1):
                    for j in range(phi.shape[0]):

                        # Note: this construction is here, because the entries of phi that belong to other DBCs should not be set to 0
                        # This specifically happens if DBCs are applied along an edge.
                        if j == cdof:
                            assert np.isclose(phi[j,:i], 0).all()
                            assert np.isclose(phi[j,i+1:], 0).all()

                            phi[j,i] = 1.0

                        elif not j in cdofs:
                            phi[j,i] = 0.0

        self._Phic = phi
        globdat['Phic'] = self._Phic

        # Get the observed force vector
        self._y = self._Phic.T @ (self._fc - self._Kc @ self._m)

        # Set alpha to the optimal value if necessary
        if self._alpha == 'opt':
            self._alpha = self._get_alpha_opt()
            print('Setting alpha to optimal value: {:.4f}'.format(self._alpha))


    def _get_prior_mean(self, params, globdat):

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')

        if field == 'u':

            ###################
            # PRIOR MEAN ON u #
            ###################

            # Return the prior of the displacement field
            params[gppn.PRIORMEAN] = self._m

        elif field == 'f':

            ###################
            # PRIOR MEAN ON f #
            ###################

            # Return the prior of the force field
            params[gppn.PRIORMEAN] = self._Kc @ self._m

        else:

            raise ValueError(field)


    def _get_posterior_mean(self, params, globdat):

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._alpha**2 * self._Phic.T @ self._Mc @ self._Phic + np.identity(self._nobs) * self._noise2

        # Do the cholesky decomposition of Sigma_obs only if necessary
        if not '_sqrtObs' in vars(self):
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

        # Solve the system for the observed forces
        if not '_v0' in vars(self):
            self._v0 = np.linalg.solve(self._sqrtObs, self._y)

        # Get the posterior of the force field
        if not '_f_post' in vars(self):
            self._f_post = self._Kc @ self._m + self._alpha**2 * self._Mc @ self._Phic @ np.linalg.solve(self._sqrtObs.T, self._v0)

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
                var_prior = np.zeros(self._dc)
                for i in range(self._dc):
                    var_prior[i] = self._alpha**2 * self._V2[i,:] @ self._V2[i,:].T

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
                var_prior = self._alpha**2 * self._Mc.diagonal()

        else:

            raise ValueError(field)

        # Return either the full covariance matrix, or only its diagonal
        if fullSigma:
            params[gppn.PRIORCOVARIANCE] = Sigma_prior
        else:
            params[gppn.PRIORCOVARIANCE] = var_prior


    def _get_posterior_covariance(self, params, globdat):

        # Check if the posterior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')
        fullSigma = params.get(gppn.FULLCOVARIANCE, False)

        # Check if the prior variance is given in the params, otherwise, take an action to obtain it
        if not gppn.PRIORCOVARIANCE in params:
            self.take_action(gpact.GETPRIORCOVARIANCE, params, globdat)

        # Get the prior variance from the params
        var_prior = params[gppn.PRIORCOVARIANCE]

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._alpha**2 * self._Phic.T @ self._Mc @ self._Phic + np.identity(self._nobs) * self._noise2

        # Do the cholesky decomposition of Sigma_obs only if necessary
        if not '_sqrtObs' in vars(self):
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

        # Solve the inverse observation covariance once for each observation
        if not '_V1' in vars(self):
            self._V1 = np.linalg.solve(self._sqrtObs, self._Phic.T @ self._Mc)

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
                Sigma_post = var_prior.copy()
                Sigma_post -= self._alpha**4 * self._V3 @ self._V3.T

            else:

                # If not, compute only the diagonal of the covariance matrix
                var_post = var_prior.copy()
                for i in range(self._dc):
                    var_post[i] -= self._alpha**4 * self._V3[i,:] @ self._V3[i,:].T

        elif field == 'f':

            #############################
            # POSTERIOR COVARIANCE ON f #
            #############################

            # Sigma = alpha2 * M - alpha4 * M * phi * inv(alpha2 * phi.T * M * phi + Sigma_e) * phi.T * M

            # Check if the full covariance matrix should be returned
            if fullSigma:

                # If so, compute the full covariance matrix
                Sigma_post = var_prior.copy()
                Sigma_post -= self._alpha**4 * self._V1.T @ self._V1

            else:

                # If not, compute only the diagonal of the covariance matrix
                var_post = var_prior.copy()
                for i in range(self._dc):
                    var_post[i] -= self._alpha**4 * self._V1[:,i].T @ self._V1[:,i]

        else:

            raise ValueError(field)

        # Return either the full covariance matrix, or only its diagonal
        if fullSigma:
            params[gppn.POSTERIORCOVARIANCE] = Sigma_post
        else:
            params[gppn.POSTERIORCOVARIANCE] = var_post


    def _get_prior_samples(self, params, globdat):

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')
        nsamples = params.get(gppn.NSAMPLE, 1)
        rng = params.get(gppn.RNG, np.random.default_rng())

        # Do the cholesky decomposition of M only if necessary
        if not '_sqrtM' in vars(self):
            self._sqrtM = np.linalg.cholesky(self._Mc)

        # Define the array that stores the samples
        samples = np.zeros((self._dc, nsamples))

        # Get a loop to create all samples
        for i in range(nsamples):

            # Get a sample from the standard normal distribution
            z = rng.standard_normal(self._dc)

            # Get the sample of the force field
            f = self._alpha * self._sqrtM @ z

            if field == 'u':

                #####################
                # PRIOR SAMPLE ON u #
                #####################

                # Compute the corresponding displacement field
                u = np.linalg.solve(self._Kc, f) + self._m

                samples[:,i] = u

            elif field == 'f':

                #####################
                # PRIOR SAMPLE ON f #
                #####################

                samples[:,i] = f + self._Kc @ self._m

            else:

                raise ValueError(field)

        # Return the array of samples
        params[gppn.PRIORSAMPLES] = samples


    def _get_posterior_samples(self, params, globdat):

        # Check if the posterior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')
        nsamples = params.get(gppn.NSAMPLE, 1)
        rng = params.get(gppn.RNG, np.random.default_rng())

        # Do the cholesky decomposition of M only if necessary
        if not '_sqrtM' in vars(self):
            self._sqrtM = np.linalg.cholesky(self._Mc)

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._alpha**2 * self._Phic.T @ self._Mc @ self._Phic + np.identity(self._nobs) * self._noise2

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
            x1 = self._alpha * self._sqrtM @ z1
            x2 = self._sqrtNoise @ z2

            # Compute the perturbed observation
            f_pert = self._y - self._Phic.T @ x1 + x2

            # Multiply the perturbed observation by the Kalman gain
            f = np.linalg.solve(self._sqrtObs, f_pert)
            f = np.linalg.solve(self._sqrtObs.T, f)
            f = self._alpha**2 * self._Mc @ self._Phic @ f

            # Add the prior sample
            f += x1

            if field == 'u':

                #########################
                # POSTERIOR SAMPLE ON u #
                #########################

                # Compute the corresponding displacement field
                u = np.linalg.solve(self._Kc, f) + self._m

                samples[:,i] = u

            elif field == 'f':

                #########################
                # POSTERIOR SAMPLE ON f #
                #########################

                samples[:,i] = f + self._Kc @ self._m

            else:

                raise ValueError(field)

        # Return the array of samples
        params[gppn.POSTERIORSAMPLES] = samples


    def _get_alpha_opt(self):

        #################
        # OPTIMAL ALPHA #
        #################

        # If so, determine the optimal value of alpha
        L = np.linalg.cholesky(self._Phic.T @ self._Mc @ self._Phic)
        v = np.linalg.solve(L, self._y)
        alpha2 = v.T @ v / self._nobs

        return np.sqrt(alpha2)


    def _get_log_likelihood(self, params, globdat):

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._alpha**2 * self._Phic.T @ self._Mc @ self._Phic + np.identity(self._nobs) * self._noise2

        # Do the cholesky decomposition of Sigma_obs only if necessary
        if not '_sqrtObs' in vars(self):
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

        # Solve the system for the observed forces
        if not '_v0' in vars(self):
            self._v0 = np.linalg.solve(self._sqrtObs, self._y)

        ##################
        # LOG LIKELIHOOD #
        ##################

        # l = - 1/2 * y * inv(Sigma + Sigma_e) * y - 1/2 * log|Sigma + Sigma_e| - n/2 * log(2*pi)

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

        # if 'dx' not in table:
        #     table['dx'] = np.zeros(len(globdat[gn.NSET]))
        # if len(self._dof_types) > 1:
        #     if 'dy' not in table:
        #         table['dy'] = np.zeros(len(globdat[gn.NSET]))
        # if len(self._dof_types) > 2:
        #     if 'dz' not in table:
        #         table['dz'] = np.zeros(len(globdat[gn.NSET]))

        # Define a dictionary for the settings of u
        u_params = {}
        u_params[gppn.FIELD] = 'u'
        u_params[gppn.FULLCOVARIANCE] = False

        # Define a dictionary for the settings of f
        f_params = {}
        f_params[gppn.FIELD] = 'f'
        f_params[gppn.FULLCOVARIANCE] = False

        if name == 'var_f_prior':
            self._get_prior_covariance(f_params, globdat)
            var = f_params[gppn.PRIORCOVARIANCE]
        elif name == 'var_f_post':
            self._get_posterior_covariance(f_params, globdat)
            var = f_params[gppn.POSTERIORCOVARIANCE]
        elif name == 'var_u_prior':
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

                # table['dx'][inode] = var[idofs[0]]
                # if len(self._dof_types) > 1:
                #     table['dy'][inode] = var[idofs[1]]
                # if len(self._dof_types) > 2:
                #     table['dz'][inode] = var[idofs[2]]

                if any(idof in cdofs for idof in idofs):
                    tbwts[inode] = np.inf
                else:
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
    factory.declare_model('GP', GPModel)
