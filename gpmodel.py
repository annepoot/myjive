import numpy as np
import scipy.linalg as spla
import scipy.sparse.linalg as spspla

from jive.fem.names import Actions as act
from jive.fem.names import GPActions as gpact
from jive.fem.names import ParamNames as pn
from jive.fem.names import GPParamNames as gppn
from jive.fem.names import GlobNames as gn
from jive.fem.names import PropNames as prn
from jive.model.model import Model
from jive.solver.constrainer import Constrainer

OBSNOISE = 'obsNoise'
PDNOISE = 'pdNoise'
PRIOR = 'prior'
TYPE = 'type'
FUNC = 'func'
PREPROJECT = 'preproject'
HYPERPARAMS = 'hyperparams'
SHAPE = 'shape'
INTSCHEME = 'intScheme'


class GPModel(Model):

    def take_action(self, action, params, globdat):
        showmsg = True

        if action == gpact.CONFIGUREFEM:
            self._configure_fem(params, globdat)
        elif action == gpact.CONFIGUREPRIOR:
            self._configure_prior(params, globdat)
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
        elif action == gpact.KALMANUPDATE:
            self._kalman_update(params, globdat)
        elif action == gpact.GETLOGLIKELIHOOD:
            self._get_log_likelihood(params, globdat)
        elif action == gpact.PROJECTSAMPLES:
            self._project_samples(params, globdat)
        elif action == act.GETTABLE:
            if 'var' in params[pn.TABLENAME]:
                self._get_variances(params, globdat)
            elif 'std' in params[pn.TABLENAME]:
                self._get_standard_deviations(params, globdat)
            else:
                showmsg = False
        else:
            showmsg = False

        if showmsg:
            print('GPModel taking action', action)

    def configure(self, props, globdat):
        self._dc = globdat[gn.DOFSPACE].dof_count()
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][prn.TYPE], props[SHAPE][INTSCHEME])
        self._rank = self._shape.global_rank()

        # Get the observational noise, and pd noise
        self._noise2 = float(props.get(OBSNOISE))**2
        self._pdnoise2 = float(props.get(PDNOISE, 1e-8))**2

        # Get the prior properties
        priorprops = props[PRIOR]
        self._prior = priorprops[TYPE]
        if self._prior == 'kernel':
            self._kernel = priorprops[FUNC]
        elif self._prior == 'SPDE':
            self._covariance = priorprops[FUNC]
        else:
            raise ValueError('prior has to be "kernel" or "SPDE"')

        self._hyperparams = {}
        for key, value in priorprops[HYPERPARAMS].items():
            if value != 'opt':
                self._hyperparams[key] = float(value)
            else:
                self._hyperparams[key] = 'opt'

        # Get the preprojection props
        self._preproject = bool(eval(props.get(PREPROJECT, 'False')))

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
        u = globdat.get(gn.STATE0)
        f = globdat.get(gn.EXTFORCE)
        c = globdat.get(gn.CONSTRAINTS)

        conmanK = Constrainer(c, K)
        conmanM = Constrainer(c, M)

        # Get the actual constrained stiffness matrix and force vector
        self._Mc = conmanM.get_output_matrix()
        self._Kc = conmanK.get_output_matrix()
        fc = conmanK.get_rhs(f)

        # Get the prior mean
        self._m = params[gppn.PRIORMEAN]
        self._cdofs, self._cvals = c.get_constraints()

        # Get the phi matrix
        self._Phi = self._get_phi(globdat)

        # Constrain the phi matrix based on Dirichlet BCs
        self._Phic = self._constrain_phi(self._Phi, self._cdofs)

        # Store Phi and Phic in globdat
        globdat['Phi'] = self._Phi
        globdat['Phic'] = self._Phic

        # Get the observation operator
        self._H = self._Phic.T @ self._Kc

        # Get the observed force vector
        self._g = self._Phic.T @ fc

        # Get the centered observation vector (as a deviation from the prior mean)
        self._y = self._g - self._H @ self._m

        # Get the coarse solution vector
        self._uc = self._Phic @ np.linalg.solve(self._H @ self._Phic, self._g)

    def _configure_prior(self, params, globdat):

        # Define a dictionary with relevant functions
        eval_dict = self._get_eval_dict()

        # Check if we have a kernel covariance
        if self._prior == 'kernel':

            # Get the covariance matrix by looping over all dofs
            self._Sigma = np.zeros((self._dc, self._dc))

            nodes = globdat[gn.NSET]
            dofspace = globdat[gn.DOFSPACE]

            for i in range(len(nodes)):
                icoords = nodes[i].get_coords()
                idofs = dofspace.get_dofs([i], dofspace.get_types())
                eval_dict['x0'] = icoords

                for j in range(len(nodes)):
                    jcoords = globdat[gn.NSET][j].get_coords()
                    jdofs = dofspace.get_dofs([j], dofspace.get_types())
                    eval_dict['x1'] = jcoords

                    self._Sigma[np.ix_(idofs, jdofs)] = eval(self._kernel, eval_dict)

        else:

            # Get the covariance matrix by 1 matrix evaluation
            self._Sigma = eval(self._covariance, eval_dict)

        # Apply boundary conditions to the prior
        self._m, self._Sigma = self._apply_covariance_bcs(self._m, self._Sigma)

    def _get_prior_mean(self, params, globdat):

        # Return the prior of the force field
        params[gppn.PRIORMEAN] = self._m

    def _get_posterior_mean(self, params, globdat):

        ##################
        # POSTERIOR MEAN #
        ##################

        # u_bar = Sigma * H.T * inv(H * Sigma * H.T + Sigma_e) * y

        # Get the relevant matrices
        self._get_sqrtObs()
        self._get_v0()

        # Get the posterior of the force field
        if not '_u_post' in vars(self):
            self._u_post = self._m + self._premul_Sigma(self._H.T @ self._solve_triangular(self._sqrtObs.T, self._v0, lower=False))

        # Return the posterior of the displacement field
        params[gppn.POSTERIORMEAN] = self._u_post

    def _get_prior_covariance(self, params, globdat):

        ####################
        # PRIOR COVARIANCE #
        ####################

        # Compute the full covariance matrix
        Sigma_prior = self._Sigma
        params[gppn.PRIORCOVARIANCE] = Sigma_prior

    def _get_posterior_covariance(self, params, globdat):

        ########################
        # POSTERIOR COVARIANCE #
        ########################

        # Sigma_bar = Sigma - Sigma * H.T * inv(H * Sigma * H.T + Sigma_e) * H * Sigma

        # Check if the prior variance is given in the params, otherwise, take an action to obtain it
        if not gppn.PRIORCOVARIANCE in params:
            self.take_action(gpact.GETPRIORCOVARIANCE, params, globdat)

        # Get the prior variance from the params
        var_prior = params[gppn.PRIORCOVARIANCE]

        # Get the relevant matrices
        self._get_V1()

        # Compute the full covariance matrix
        Sigma_post = var_prior.copy()
        Sigma_post -= self._V1.T @ self._V1
        params[gppn.POSTERIORCOVARIANCE] = Sigma_post

    def _get_prior_samples(self, params, globdat):

        nsamples = params.get(gppn.NSAMPLE, 1)
        rng = params.get(gppn.RNG, np.random.default_rng())

        #################
        # PRIOR SAMPLES #
        #################

        # u = m + sqrt(Sigma) * z

        # Get the relevant matrices
        self._get_sqrtSigma()

        # Define the array that stores the samples
        samples = np.zeros((self._dc, nsamples))

        # Get a loop to create all samples
        for i in range(nsamples):

            # Get a sample from the standard normal distribution
            # NOTE: self._sqrtSigma.shape[1] does not have to be the same as self._dc
            # For example, if Ensemble Kalman is used, it is equal to self._nens instead.
            z = rng.standard_normal(self._sqrtSigma.shape[1])

            # Perform an initial projection over the coarse space
            if self._preproject:
                self._get_P()
                z = self._P @ z

            # Get the sample of the force field
            u = self._sqrtSigma @ z + self._m

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

        # Get the relevant matrices
        self._get_sqrtObs()
        self._get_sqrtSigma()
        self._get_sqrtNoise()

        # Define the array that stores the samples
        samples = np.zeros((self._dc, nsamples))

        # Get a loop to create all samples
        for i in range(nsamples):

            # Get two samples from the standard normal distribution
            # Get a sample from the standard normal distribution
            # NOTE: self._sqrtSigma.shape[1] does not have to be the same as self._dc
            # For example, if Ensemble Kalman is used, it is equal to self._nens instead.
            z1 = rng.standard_normal(self._sqrtSigma.shape[1])
            z2 = rng.standard_normal(self._nobs)

            # Perform an initial projection over the coarse space
            if self._preproject:
                self._get_P()
                z1 = self._P @ z1

            # Get x1 ~ N(0, alpha^2 M) and x2 ~ N(0, Sigma_e)
            x1 = self._sqrtSigma @ z1
            x2 = self._sqrtNoise @ z2

            # Compute the perturbed observation
            u_pert = self._y - self._H @ x1 + x2

            # Multiply the perturbed observation by the Kalman gain
            u = self._solve_triangular(self._sqrtObs, u_pert, lower=True)
            u = self._solve_triangular(self._sqrtObs.T, u, lower=False)
            u = self._premul_Sigma(self._H.T @ u)

            # Add the prior sample and mean
            u += x1 + self._m

            samples[:,i] = u

        # Return the array of samples
        params[gppn.POSTERIORSAMPLES] = samples

    def _kalman_update(self, params, globdat):

        samples_prior = params[gppn.PRIORSAMPLES]
        rng = params.get(gppn.RNG, np.random.default_rng())

        #####################
        # POSTERIOR SAMPLES #
        #####################

        # u = m + sqrt(Sigma) * z1 + Sigma * H.T * inv(H * Sigma * H + Sigma_e) * (y - H * sqrt(Sigma) * z1 + sqrt(Sigma_e) * z2)

        # Get the relevant matrices
        self._get_sqrtObs()
        self._get_sqrtNoise()

        # Define the array that stores the samples
        samples_post = np.zeros(samples_prior.shape, dtype=samples_prior.dtype)

        # Get a loop to create all samples
        for i, u_prior in enumerate(samples_prior.T):

            # Center the prior sample
            u_centered = u_prior - self._m

            # Get x ~ N(0, Sigma_e)
            z = rng.standard_normal(self._nobs)
            x = self._sqrtNoise @ z

            # Compute the perturbed observation
            u_pert = self._y - self._H @ u_centered + x

            # Multiply the perturbed observation by the Kalman gain
            tmp = self._solve_triangular(self._sqrtObs, u_pert, lower=True)
            tmp = self._solve_triangular(self._sqrtObs.T, tmp, lower=False)
            u_post = self._premul_Sigma(self._H.T @ tmp)

            # Add the prior sample and mean
            u_post += u_centered + self._m

            samples_post[:,i] = u_post

        # Return the array of samples
        params[gppn.POSTERIORSAMPLES] = samples_post

    def _project_samples(self, params, globdat):

        # Get the projection matrix
        self._get_P()

        # Project the prior samples
        params[gppn.PRIORSAMPLES] = self._P @ params[gppn.PRIORSAMPLES]

        # Project the posterior samples
        params[gppn.POSTERIORSAMPLES] = self._P @ params[gppn.POSTERIORSAMPLES]

    def _get_log_likelihood(self, params, globdat):

        ##################
        # LOG LIKELIHOOD #
        ##################

        # l = - 1/2 * y * inv(H * Sigma * H.T + Sigma_e) * y - 1/2 * log|H * Sigma * H.T + Sigma_e| - n/2 * log(2*pi)

        # Get the relevant matrices
        self._get_sqrtObs()
        self._get_v0()

        l = - 0.5 * self._v0.T @ self._v0 - np.sum(np.log(self._sqrtObs.diagonal())) - 0.5 * self._nobs * np.log(2*np.pi)

        params[gppn.LOGLIKELIHOOD] = l

    def _get_eval_dict(self):

        # Define a dictionary with relevant functions
        eval_dict = {'inv':np.linalg.inv, 'exp':np.exp, 'norm':np.linalg.norm, 'np':np}
        eval_dict.update(self._hyperparams)

        # Check if we have an SPDE covariance
        if self._prior == 'SPDE':

            # Add the mass and stiffness matrices to the dictionary
            eval_dict['M'] = self._Mc.toarray()
            eval_dict['K'] = self._Kc.toarray()
            eval_dict['Phi'] = self._Phi

        return eval_dict

    def _get_phi(self, globdat):

        elems = globdat[gn.ESET]
        elemsc = globdat[gn.COARSEMESH][gn.ESET]
        nodes = globdat[gn.NSET]
        nodesc = globdat[gn.COARSEMESH][gn.NSET]
        dofs = globdat[gn.DOFSPACE]
        dofsc = globdat[gn.COARSEMESH][gn.DOFSPACE]

        Phi = np.zeros((dofs.dof_count(), dofsc.dof_count()))

        # Go over the coarse mesh
        for elemc in elemsc:
            inodesc = elemc.get_nodes()
            idofsc = dofsc.get_dofs(inodesc, self._dof_types)
            coordsc = nodesc.get_some_coords(inodesc)

            # Get the bounding box of the coarse element
            bbox = np.zeros((self._rank, 2))
            for i in range(self._rank):
                bbox[i,0] = min(coordsc[i,:])
                bbox[i,1] = max(coordsc[i,:])

            # Go over the fine mesh
            for elem in elems:
                inodes = elem.get_nodes()
                coords = nodes.get_some_coords(inodes)

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
                            svals = np.round(self._shape.eval_shape_functions(loc_point), 12)

                            # Get the dofs belonging to this node
                            idofs = dofs.get_dofs([inodes[n]], self._dof_types)

                            # Store the relative shape function values in the phi matrix
                            for i, idof in enumerate(idofs):
                                Phi[idof, idofsc[i::len(self._dof_types)]] = svals

        return Phi

    def _constrain_phi(self, Phi, cdofs):
        Phic = Phi.copy()

        for i in range(Phic.shape[1]):
            for cdof in cdofs:
                if np.isclose(Phic[cdof,i], 1):
                    Phic[:,i] = 0.0
                    Phic[cdof,i] = 1.0

        return Phic

    def _apply_covariance_bcs(self, m, Sigma):
        Sigmac = Sigma.copy()
        mc = m.copy()

        # Split Sigma along boundary and internal nodes
        idofs = np.delete(np.arange(self._dc), self._cdofs)

        Sigma_bb = Sigma[np.ix_(self._cdofs,self._cdofs)]
        Sigma_bi = Sigma[np.ix_(self._cdofs,idofs)]
        Sigma_ib = Sigma[np.ix_(idofs,self._cdofs)]

        Sigma_bb_inv = np.linalg.inv(Sigma_bb)

        # Update the prior mean by observing the displacement at the bcs
        mc[idofs] += Sigma_ib @ Sigma_bb_inv @ (self._cvals - m[self._cdofs])
        mc[self._cdofs] = self._cvals

        # Update the prior covariance as well
        Sigmac[np.ix_(idofs,idofs)] -= Sigma_ib @ Sigma_bb_inv @ Sigma_bi

        # Decouple the bc covariance from the internal nodes
        Sigmac[self._cdofs,:] = Sigmac[:,self._cdofs] = 0.0

        # Add a tiny noise to ensure Sigma is positive definite rather than semidefinite
        Sigmac += self._pdnoise2 * np.identity(self._dc)

        return mc, Sigmac

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

    def _solve_triangular(self, A, b, lower):
        return spla.solve_triangular(A, b, lower=lower)

    def _get_Sigma_obs(self):

        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._H @ self._Sigma @ self._H.T + np.identity(self._nobs) * self._noise2

    def _get_sqrtObs(self):

        if not '_sqrtObs' in vars(self):
            self._get_Sigma_obs()
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

    def _get_sqrtSigma(self):

        if not '_sqrtSigma' in vars(self):
            self._sqrtSigma = np.linalg.cholesky(self._Sigma)

    def _get_sqrtNoise(self):

        if not '_sqrtNoise' in vars(self):
            self._sqrtNoise = np.sqrt(self._noise2) * np.identity(self._nobs)

    def _get_v0(self):

        if not '_v0' in vars(self):
            self._get_sqrtObs()
            self._v0 = self._solve_triangular(self._sqrtObs, self._y, lower=True)

    def _get_V1(self):

        if not '_V1' in vars(self):
            self._get_sqrtObs()
            self._V1 = self._postmul_Sigma(self._solve_triangular(self._sqrtObs, self._H, lower=True))

    def _get_P(self):

        if not '_P' in vars(self):
            self._P = self._Phi @ np.linalg.inv(self._Phi.T @ self._Phi) @ self._Phi.T

    def _premul_Sigma(self, X):

        return self._Sigma @ X

    def _postmul_Sigma(self, X):

        return X @ self._Sigma


def declare(factory):
    factory.declare_model('GP', GPModel)
