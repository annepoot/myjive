__all__ = ["GPActions", "GPParamNames"]


class GPActions:
    CONFIGUREFEM = "configureFEM"
    CONFIGUREPRIOR = "configurePrior"
    GETPRIORMEAN = "getPriorMean"
    GETPRIORCOVARIANCE = "getPriorCovariance"
    GETPOSTERIORMEAN = "getPosteriorMean"
    GETPOSTERIORCOVARIANCE = "getPosteriorCovariance"
    GETPRIORSAMPLES = "getPriorSamples"
    GETPOSTERIORSAMPLES = "getPosteriorSamples"
    GETLOGLIKELIHOOD = "getLogLikelihood"
    KALMANUPDATE = "kalmanUpdate"
    PROJECTSAMPLES = "projectSamples"


class GPParamNames:
    FIELD = "field"
    PRIORMEAN = "priorMean"
    PRIORCOVARIANCE = "priorCovariance"
    POSTERIORMEAN = "posteriorMean"
    POSTERIORCOVARIANCE = "posteriorCovariance"
    PRIORSAMPLES = "priorSamples"
    POSTERIORSAMPLES = "posteriorSamples"
    LOGLIKELIHOOD = "logLikelihood"
    FULLCOVARIANCE = "fullCovariance"
    NSAMPLE = "nsamples"
    SEED = "seed"
    RNG = "randomNumberGenerator"
