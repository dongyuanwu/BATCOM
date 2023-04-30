
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_psi.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::export]]
arma::vec mvrnormArma_mean0(arma::mat const &sigma) {
    int ncols = sigma.n_cols;
    arma::vec Y = arma::randn(ncols);
    return arma::chol(sigma, "lower") * Y;
}

// [[Rcpp::export]]
arma::ivec oneMultinomCalt(arma::vec &probs) {
    int k = probs.size();
    arma::ivec ans(k);
    arma::vec probs_transform = probs / arma::accu(probs);
    rmultinom(1, probs_transform.begin(), k, ans.begin());
    return(ans);
}

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::vec &x, arma::mat const &trimat){
    arma::uword const n = trimat.n_cols;

    for(unsigned j = n; j-- > 0;){
        double tmp(0.);
        for(unsigned i = 0; i <= j; ++i)
            tmp += trimat.at(i, j) * x[i];
        x[j] = tmp;
    }
}

static double const log2pi = std::log(2.0 * M_PI);


// [[Rcpp::export]]
double dmvnrm_arma_fast_vec(arma::vec const &x,  
                             arma::vec const &sigma, 
                             bool const logd = false) { 
    using arma::uword;
    uword const xdim = x.n_elem;
    double const rootisum = - arma::accu(log(sigma)) / 2.0, 
        constants = -(double)xdim/2.0 * log2pi, 
        other_terms = rootisum + constants;
    
    double out = other_terms - 0.5 * arma::accu(x % x / sigma);    
    
    if (logd)
        return out;
    return exp(out);
}


// [[Rcpp::export]]
double dmvnrm_diff(arma::vec const &x1, arma::vec const &x2, arma::mat const &sigma, 
                   bool const logd = false) { 
    using arma::uword;
    double out;
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    
    arma::vec z1 = x1, z2 = x2;
    inplace_tri_mat_mult(z1, rooti);
    inplace_tri_mat_mult(z2, rooti);
    out = - 0.5 * (arma::dot(z1, z1) - arma::dot(z2, z2));    
    
    if (logd)
        return out;
    return exp(out);
}


// [[Rcpp::export]]
arma::mat logdens_ddparams_re_sp(arma::vec const &params, arma::sp_mat const &X_beta_re, arma::vec const &Y,
                           arma::vec const &eta, arma::vec const &TT, int const &nspot, arma::vec const &beta_re_cov) {
    // Transform and check parameters
    
    int q_beta_re = X_beta_re.n_cols;
    int q = q_beta_re - 2*nspot;
    arma::sp_mat X(X_beta_re.cols(0, q-1));
    
    double logphi = params[0];
    double theta = params[1];
    if (theta > 10) {
        theta = 10;
    } else if (theta < -10) {
        theta = -10;
    }
    double exptheta = exp(theta);
    double p = (2*exptheta + 1) / (exptheta + 1);
    arma::vec beta_re = params.rows(2, params.n_elem-1);
    double p1 = p - 1;
    double p2 = 2 - p;
    
    arma::uvec ids = find(Y > 0); // Find indices
    arma::vec subTT = TT.rows(ids);
    arma::vec subY = Y.rows(ids);
    
    arma::vec subeta = eta.rows(ids);
    
    arma::sp_mat dd_bb(q_beta_re, q_beta_re);
    arma::vec tempval2 = exp(p2*eta - logphi);
    arma::vec tempval1 = exp(-p1*eta - logphi);
    arma::vec subtempval1 = tempval1.rows(ids);
    dd_bb = - (X_beta_re.t() * arma::diagmat(Y % tempval1 * p1 + tempval2 * p2) * X_beta_re);
    dd_bb.diag() = dd_bb.diag() - 1 / beta_re_cov;
    
    double dd_logphi;
    double alpha = p2 / p1;
    dd_logphi = - arma::accu(tempval2 / p2 + Y % tempval1 / p1);
    dd_logphi = dd_logphi - 1 / 100;
    double p_dtheta = p1 * p2;
    double dd_theta;
    
    dd_theta = arma::accu(tempval2 % eta * p2 - tempval2 - eta % eta % tempval2 * p2 * p_dtheta) * exptheta;
    arma::vec psivec(subTT.size(), arma::fill::zeros);
    arma::vec dpsivec(subTT.size(), arma::fill::zeros);
    
    for (unsigned int i = 0; i < subTT.size(); i++) {
        psivec[i] = gsl_sf_psi(subTT[i]*alpha);
        dpsivec[i] = gsl_sf_psi_1(subTT[i]*alpha);
    }
    
    dd_theta = dd_theta - arma::accu(subY % (subeta % subeta % subtempval1 * p1 * p_dtheta + 
        subeta % subtempval1 * p1 + subtempval1)) / exptheta;
    dd_theta = dd_theta + arma::accu(subTT * (p_dtheta + alpha * p2 + p2 * p2 - alpha * logphi + alpha * p2) -
        alpha * subTT % (psivec + log(p1) - log(subY) + alpha * subTT % dpsivec));
    dd_theta = dd_theta - 4 * exptheta / (1 + exptheta) / (1 + exptheta);
    
    double dd_logphi_theta;
    arma::vec dd_logphi_bb(q, arma::fill::zeros), dd_theta_bb(q, arma::fill::zeros);
    
    
    dd_logphi_theta = exptheta * arma::accu(tempval2 - tempval2 % eta * p2);
    dd_logphi_theta = dd_logphi_theta - arma::accu(subY % (subtempval1 % subeta * p1 + subtempval1)) / exptheta + arma::accu(subTT) * alpha;
    
    dd_logphi_bb = X_beta_re.t() * (tempval2 - Y % tempval1);
    dd_theta_bb = X_beta_re.t() * ((tempval2 - Y % tempval1) % eta) * p_dtheta;
    
    
    arma::sp_mat ddparams(q_beta_re+2, q_beta_re+2);
    ddparams(0, 0) = dd_logphi;
    ddparams(1, 1) = dd_theta;
    ddparams(0, 1) = dd_logphi_theta;
    ddparams(1, 0) = dd_logphi_theta;
    ddparams.submat(0, 2, 0, q_beta_re+1) = dd_logphi_bb.t();
    ddparams.submat(1, 2, 1, q_beta_re+1) = dd_theta_bb.t();
    ddparams.submat(2, 0, q_beta_re+1, 0) = dd_logphi_bb;
    ddparams.submat(2, 1, q_beta_re+1, 1) = dd_theta_bb;
    ddparams.submat(2, 2, q_beta_re+1, q_beta_re+1) = dd_bb;
    
    return (arma::mat(ddparams));
    
}


// [[Rcpp::export]]
arma::mat logdens_ddbeta(arma::vec const &beta, arma::mat const &X, 
                         arma::vec const &Y, arma::vec const &beta_cov, 
                         double const &logphi, double const &p) {
    
    // Transform and check parameters
    
    int q = beta.size();
    
    double p2 = 2 - p;
    double p1 = p - 1;
    
    arma::uvec ids = find(Y > 0); // Find indices
    arma::vec subY = Y.rows(ids);
    
    arma::vec eta = X * beta;
    arma::vec subeta = eta.rows(ids);
    
    arma::mat dd_bb(q, q, arma::fill::zeros);
    arma::vec tempval2 = exp(p2*eta - logphi);
    arma::vec tempval1 = exp(-p1*eta - logphi);
    arma::vec subtempval1 = tempval1.rows(ids);
    dd_bb = - (X.t() * arma::diagmat(Y % tempval1 * p1 + tempval2 * p2) * X);
    dd_bb.diag() = dd_bb.diag() - 1 / beta_cov;
    
    return (dd_bb);
    
}

// [[Rcpp::export]]
arma::sp_mat logdens_ddbeta_re_sp(arma::vec const &beta_re, arma::sp_mat const &X_beta_re, 
                                  arma::vec const &Y, arma::vec const &tempval1, 
                                  arma::vec const &tempval2, double const &p, 
                                  arma::vec const &beta_re_cov) {
    
    // Transform and check parameters
    
    int q_beta_re = beta_re.size();
    
    double p2 = 2 - p;
    double p1 = p - 1;
    
    arma::sp_mat dd_bb(q_beta_re, q_beta_re);
    
    dd_bb = - (X_beta_re.t() * arma::diagmat(Y % tempval1 * p1 + tempval2 * p2) * X_beta_re);
    dd_bb.diag() = dd_bb.diag() - 1 / beta_re_cov;
    
    return (dd_bb);
    
}


// [[Rcpp::export]]
arma::vec logdens_dparams(arma::vec const &params, arma::mat const &X, arma::vec const &Y, 
                          arma::vec const &TT, arma::vec const &beta_cov,
                          arma::vec const &reL, arma::vec const &reR) {
    
    // Transform and check parameters
    
    int q = params.size() - 2;
    
    double logphi = params[0];
    double theta = params[1];
    if (theta > 10) {
        theta = 10;
    } else if (theta < -10) {
        theta = -10;
    }
    double exptheta = exp(theta);
    double p = (2*exptheta + 1) / (exptheta + 1);
    arma::vec beta = params.rows(2, q+1);
    double p2 = 2 - p;
    double p1 = p - 1;
    
    arma::uvec ids = find(Y > 0); // Find indices
    arma::vec subTT = TT.rows(ids);
    arma::vec subY = Y.rows(ids);
    
    arma::vec eta = X * beta + reL + reR;
    for (unsigned int j = 0; j < eta.n_rows; j++) {
        if (eta[j] > 10) eta[j] = 10;
    }
    arma::vec subeta = eta.rows(ids);
    
    arma::vec d_bb(q, arma::fill::zeros);
    arma::vec tempval2 = exp(p2*eta - logphi);
    arma::vec tempval1 = exp(-p1*eta - logphi);
    d_bb = X.t() * (Y % tempval1 - tempval2);
    d_bb = d_bb - beta / beta_cov;
    
    double d_logphi;
    double alpha = p2 / p1;
    d_logphi = arma::accu(tempval2 / p2 - TT - TT*alpha + Y % tempval1 / p1) + 1;
    d_logphi = d_logphi - logphi / 100;
    double d_theta;
    
    d_theta = arma::accu(tempval2 % eta * p2 - tempval2) * exptheta;
    arma::vec psivec(subTT.size(), arma::fill::zeros);
    
    for (unsigned int i = 0; i < subTT.size(); i++) {
        psivec[i] = gsl_sf_psi(subTT[i]*alpha);
    }
    d_theta = d_theta + arma::accu(subTT * p1 - subTT * alpha * p2 +
        (psivec % subTT + subTT * log(p1) - subTT % log(subY) + subTT * logphi) * alpha +
        subY / exptheta % (tempval1.rows(ids) % subeta * p1 + tempval1.rows(ids))) +
        (1 - exptheta) / (1 + exptheta);
    d_theta = d_theta + (1 - exptheta) / (1 + exptheta);
    
    arma::vec dparams(q+2, arma::fill::zeros);
    dparams[0] = d_logphi;
    dparams[1] = d_theta;
    dparams.subvec(2, q+1) = d_bb;
    
    return (dparams);
    
}

// [[Rcpp::export]]
arma::vec logdens_dparams_re_sp(arma::vec const &params, arma::sp_mat const &X_beta_re, 
                                arma::vec const &Y, arma::vec const &eta, 
                                arma::vec const &TT, int const &nspot, 
                                arma::vec const &beta_re_cov) {
    
    // Transform and check parameters
    
    int q_beta_re = X_beta_re.n_cols;
    int q = q_beta_re - 2*nspot;
    arma::sp_mat X(X_beta_re.cols(0, q-1));
    
    double logphi = params[0];
    double theta = params[1];
    if (theta > 10) {
        theta = 10;
    } else if (theta < -10) {
        theta = -10;
    }
    double exptheta = exp(theta);
    double p = (2*exptheta + 1) / (exptheta + 1);
    arma::vec beta_re = params.rows(2, params.n_elem-1);
    double p1 = p - 1;
    double p2 = 2 - p;
    
    
    arma::uvec ids = find(Y > 0); // Find indices
    arma::vec subTT = TT.rows(ids);
    arma::vec subY = Y.rows(ids);
    
    arma::vec subeta = eta.rows(ids);
    
    arma::vec d_bb(q_beta_re, arma::fill::zeros);
    arma::vec tempval2 = exp(p2*eta - logphi);
    arma::vec tempval1 = exp(-p1*eta - logphi);
    d_bb = X_beta_re.t() * (Y % tempval1 - tempval2);
    d_bb = d_bb - beta_re / beta_re_cov;
    
    double d_logphi;
    double alpha = p2 / p1;
    d_logphi = arma::accu(tempval2 / p2 - TT - TT*alpha + Y % tempval1 / p1) + 1;
    d_logphi = d_logphi - logphi / 100;
    double d_theta;
    
    d_theta = arma::accu(tempval2 % eta * p2 - tempval2) * exptheta;
    arma::vec psivec(subTT.size(), arma::fill::zeros);
    
    for (unsigned int i = 0; i < subTT.size(); i++) {
        psivec[i] = gsl_sf_psi(subTT[i]*alpha);
    }
    d_theta = d_theta + arma::accu(subTT * p1 - subTT * alpha * p2 +
        (psivec % subTT + subTT * log(p1) - subTT % log(subY) + subTT * logphi) * alpha +
        subY / exptheta % (tempval1.rows(ids) % subeta * p1 + tempval1.rows(ids))) +
        (1 - exptheta) / (1 + exptheta);
    d_theta = d_theta + (1 - exptheta) / (1 + exptheta);
    
    arma::vec dparams(q_beta_re+2, arma::fill::zeros);
    dparams[0] = d_logphi;
    dparams[1] = d_theta;
    dparams.subvec(2, q_beta_re+1) = d_bb;
    
    return (dparams);
    
}


// [[Rcpp::export]]
arma::vec logdens_dbeta_re_sp(arma::vec const &beta_re, arma::sp_mat const &X_beta_re, 
                           arma::vec const &Y, arma::vec const &tempval1, 
                           arma::vec const &tempval2, arma::vec const &beta_re_cov) {
    
    // Transform and check parameters
    
    int q_beta_re = X_beta_re.n_cols;
    arma::vec d_bb(q_beta_re, arma::fill::zeros);
    
    d_bb = X_beta_re.t() * (Y % tempval1 - tempval2);
    d_bb = d_bb - beta_re / beta_re_cov;
    
    return (d_bb);
    
}


// [[Rcpp::export]]
double logdens(arma::vec const &params, arma::mat const &X, arma::vec const &Y, 
               arma::vec const &TT, arma::vec const &beta_cov,
               arma::vec const &reL, arma::vec const &reR){

    // Transform and check parameters
    
    int q = params.size() - 2;
    
    double logphi = params[0];
    double phi = exp(logphi);
    double theta = params[1];
    if (theta > 10) {
        theta = 10;
    } else if (theta < -10) {
        theta = -10;
    }
    double exptheta = exp(theta);
    double p = (2*exptheta + 1) / (exptheta + 1);
    double p2 = 2 - p;
    double p1 = p - 1;
    
    arma::vec beta = params.rows(2, q+1);
    
    arma::vec eta = X * beta + reL + reR;
    for (unsigned int j = 0; j < eta.n_rows; j++) {
        if (eta[j] > 10) eta[j] = 10;
    }
    
    arma::vec lambda = exp(eta*p2)/(phi*p2);
    
    arma::uvec ids = find(Y > 0); // Find indices
    arma::uvec ids1 = find(Y == 0);
    arma::vec subTT = TT.rows(ids);
    arma::vec sublambda = lambda.rows(ids);
    arma::vec sublambda1 = lambda.rows(ids1);
    double logp = - arma::accu(sublambda1);
    for (unsigned int i = 0; i < subTT.size(); i++) logp = logp + log(gsl_ran_poisson_pdf(subTT[i], sublambda[i]));
    
    arma::vec subeta = eta.rows(ids);
    arma::vec subgamma = phi*p1*exp(p1*subeta);
    double alpha = p2 / p1;
    arma::vec subY = Y.rows(ids);
    for (unsigned int i = 0; i < subTT.size(); i++) logp = logp + log(gsl_ran_gamma_pdf(subY[i], alpha*subTT[i], subgamma[i]));
    
    logp = logp + logphi + log(p1*p2);
    
    logp = logp + dmvnrm_arma_fast_vec(beta, beta_cov, TRUE);
    logp = logp + log(gsl_ran_gaussian_pdf(logphi, 10));
    logp = logp + log(gsl_ran_logistic_pdf(theta, 1));
    
    return logp;

}



// [[Rcpp::export]]
arma::vec loglik_ind(arma::vec const &params, arma::mat const &X, arma::vec const &Y, 
                      arma::vec const &TT, arma::vec const &beta_cov,
                      arma::vec const &reL, arma::vec const &reR){
    
    // Transform and check parameters
    
    int q = params.size() - 2;
    
    double logphi = params[0];
    double phi = exp(logphi);
    double theta = params[1];
    if (theta > 10) {
        theta = 10;
    } else if (theta < -10) {
        theta = -10;
    }
    double exptheta = exp(theta);
    double p = (2*exptheta + 1) / (exptheta + 1);
    double p2 = 2 - p;
    double p1 = p - 1;
    
    arma::vec beta = params.rows(2, q+1);
    
    arma::vec eta = X * beta + reL + reR;
    for (unsigned int j = 0; j < eta.n_rows; j++) {
        if (eta[j] > 10) eta[j] = 10;
    }
    
    arma::vec lambda = exp(eta*p2)/(phi*p2);
    
    arma::uvec ids = find(Y > 0); // Find indices
    arma::vec subTT = TT.rows(ids);
    arma::vec sublambda = lambda.rows(ids);
    arma::vec logp = - lambda;
    
    arma::vec subeta = eta.rows(ids);
    arma::vec subgamma = phi*p1*exp(p1*subeta);
    double alpha = p2 / p1;
    arma::vec subY = Y.rows(ids);
    
    arma::vec poilp(subTT.size(), arma::fill::zeros);
    arma::vec gamlp(subTT.size(), arma::fill::zeros);
    for (unsigned int i = 0; i < subTT.size(); i++) {
        
        poilp[i] = log(gsl_ran_poisson_pdf(subTT[i], sublambda[i]));
        gamlp[i] = log(gsl_ran_gamma_pdf(subY[i], alpha*subTT[i], subgamma[i]));
        
    }
    logp.rows(ids) = poilp + gamlp;
    
    // logp = logp + logphi + log(p1*p2);
    
    // logp = logp + dmvnrm_arma_fast_vec(beta, beta_cov, TRUE);
    // logp = logp + log(gsl_ran_gaussian_pdf(logphi, 10));
    // logp = logp + log(gsl_ran_logistic_pdf(theta, 1));
    
    return logp;
    
}



// [[Rcpp::export]]
Rcpp::List newton_raphson_re_sp(arma::vec const &params_init, arma::sp_mat const &X_beta_re, arma::vec const &Y, 
                            arma::vec const &beta_re_cov, int const &nspot, arma::vec const &Lposi, 
                            arma::vec const &Rposi, double const &diff, int const &max_iter, double const &eps) {
    
    // Transform and check parameters
    
    int q_beta_re = X_beta_re.n_cols; int n = X_beta_re.n_rows;
    int q = q_beta_re - 2*nspot;
    arma::sp_mat X(X_beta_re.cols(0, q-1));
    
    arma::vec TT = arma::zeros(n);
    arma::uvec ids = find(Y > 0); // Find indices
    
    arma::vec params(q_beta_re+2, arma::fill::zeros), params_new(q_beta_re+2, arma::fill::zeros);
    params_new = params_init;
    
    arma::vec firstd(q_beta_re+2, arma::fill::zeros);
    arma::mat secondd(q_beta_re+2, q_beta_re+2, arma::fill::zeros);
    arma::sp_mat Imat = arma::speye<sp_mat>(q_beta_re+2, q_beta_re+2);
    
    int i = 0;
    
    while (arma::max(abs(params - params_new)) > diff && i < max_iter) {
        
        params = params_new;
        arma::vec beta = params.rows(2, q+1);
        arma::vec eta = X * beta;
        double logphi = params[0];
        double phi = exp(logphi);
        double theta = params[1];
        if (theta > 10) {
            theta = 10;
        } else if (theta < -10) {
            theta = -10;
        }
        double exptheta = exp(theta);
        double p = (2*exptheta + 1) / (exptheta + 1);
        
        for (int j = 0; j < nspot; j++) {
            arma::uvec Lid = find(Lposi == j);
            arma::uvec Rid = find(Rposi == j);
            eta.rows(Lid) = eta.rows(Lid) + arma::vec(Lid.n_elem, arma::fill::value(params[j+q+2]));
            eta.rows(Rid) = eta.rows(Rid) + arma::vec(Rid.n_elem, arma::fill::value(params[j+q+nspot+2]));
        }
        
        for (unsigned int j = 0; j < eta.n_rows; j++) {
            if (eta[j] > 10) eta[j] = 10;
        }
        arma::vec LAMBDA = arma::exp(eta*(2-p))/(phi*(2-p));
        arma::vec GAMMA = phi*(p-1)*exp((p-1)*eta);
        double ALPHA = (2-p)/(p-1);
        int TT_max = 5;
        arma::vec TT_prob(TT_max, arma::fill::zeros);
        for (unsigned int i = 0; i < ids.size(); i++) {
            
            for (int j = 0; j < TT_max; j++) {
                TT_prob[j] = gsl_ran_poisson_pdf((j+1), LAMBDA[ids[i]]) * gsl_ran_gamma_pdf(Y[ids[i]], ALPHA*(j+1), GAMMA[ids[i]]);
            }
            if(all(TT_prob == 0)) {
                TT[ids[i]] = TT_max;
            } else {
                TT[ids[i]] = TT_prob.index_max() + 1;
            }
        }
            
        firstd = logdens_dparams_re_sp(params, X_beta_re, Y, eta, TT, nspot, beta_re_cov);
        secondd = logdens_ddparams_re_sp(params, X_beta_re, Y, eta, TT, nspot, beta_re_cov);
        
        if (arma::det(secondd) == 0) {
            params_new = params - arma::pinv(secondd) * firstd * eps;
        } else {
            params_new = params - arma::inv(secondd, arma::inv_opts::allow_approx) * firstd * eps;
        }
        
        i = i + 1;
        
    }
    
    return Rcpp::List::create(Rcpp::Named("params_new")=params_new,
                              Rcpp::Named("TT")=TT);
    
}


// [[Rcpp::export]]
arma::vec newton_raphson_beta_re_sp(arma::sp_mat const &X_beta_re, arma::vec const &Y, 
                                    arma::vec const &beta_re_cov, double const &logphi, 
                                    double const &p, int const &nspot, arma::vec const &Lposi, 
                                    arma::vec const &Rposi, double const &diff, int const &max_iter) {
    
    // Transform and check parameters
    
    int q_beta_re = X_beta_re.n_cols;
    int q = q_beta_re - 2*nspot;
    arma::sp_mat X(X_beta_re.cols(0, q-1));
    double p1 = p - 1;
    double p2 = 2 - p;
    
    int i = 0;
    
    arma::vec beta_re(q_beta_re, arma::fill::zeros), beta_re_new(q_beta_re, arma::fill::zeros);
    beta_re_new[0] = 1;
    
    arma::vec firstd(q+2, arma::fill::zeros);
    arma::mat secondd(q+2, q+2, arma::fill::zeros);
    arma::sp_mat Imat = arma::speye<sp_mat>(q+2, q+2);
    
    while (arma::max(abs(beta_re - beta_re_new)) > diff && i < max_iter) {
        
        beta_re = beta_re_new;
        arma::vec beta = beta_re.rows(0, q-1);
        arma::vec eta = X * beta;
        
        for (int j = 0; j < nspot; j++) {
            arma::uvec Lid = find(Lposi == j);
            arma::uvec Rid = find(Rposi == j);
            eta.rows(Lid) = eta.rows(Lid) + arma::vec(Lid.n_elem, arma::fill::value(beta_re[j+q]));
            eta.rows(Rid) = eta.rows(Rid) + arma::vec(Rid.n_elem, arma::fill::value(beta_re[j+q+nspot]));
        }
        arma::vec tempval2 = exp(p2*eta - logphi);
        arma::vec tempval1 = exp(-p1*eta - logphi);
        
        firstd = logdens_dbeta_re_sp(beta_re, X_beta_re, Y, tempval1, tempval2, beta_re_cov);
        secondd = logdens_ddbeta_re_sp(beta_re, X_beta_re, Y, tempval1, tempval2, p, beta_re_cov);
        if (arma::det(secondd) == 0) {
            beta_re_new = beta_re - arma::pinv(secondd) * firstd;
        } else {
            beta_re_new = beta_re - arma::inv(secondd, arma::inv_opts::allow_approx) * firstd;
        }

        i = i + 1;
        
    }
    
    return beta_re_new;
    
}

// [[Rcpp::export]]
Rcpp::List hmc(arma::vec const &Y, arma::mat const &X, arma::vec const &init_params,
               arma::vec const &init_TT, arma::vec const &reL, arma::vec const &reR, 
               int const &niter) {
    
    gsl_rng_env_setup();
    gsl_rng *s = gsl_rng_alloc(gsl_rng_mt19937); // Create RNG seed
    gsl_rng_set(s, time(NULL)); // Seed with time
    
    int q = X.n_cols; int n = X.n_rows;
    
    arma::vec beta_cov(q, arma::fill::ones);
    beta_cov[0] = 10000;
    
    arma::vec TT = init_TT;
    arma::uvec ids = find(Y > 0); // Find indices
    
    arma::vec temp;
    
    arma::vec params = init_params;
    double logphi = params[0];
    double phi = exp(logphi);
    double theta = params[1];
    if (theta > 10) {
        theta = 10;
    } else if (theta < -10) {
        theta = -10;
    }
    double exptheta = exp(theta);
    double p = (2*exptheta + 1) / (exptheta + 1);
    arma::vec beta = params.rows(2, q+1);
    
    arma::mat Ham_M(q+2, q+2, arma::fill::eye);
    Ham_M.submat(2, 2, q+1, q+1) = - logdens_ddbeta(beta, X, Y, beta_cov, logphi, p);
    // Ham_M.submat(2, 2, q+1, q+1) = X.t() * X;
    Ham_M.submat(2, 2, q+1, q+1) = Ham_M.submat(2, 2, q+1, q+1) / Ham_M(2, 2);
    Ham_M = Ham_M * 10000;
    
    if(!Ham_M.is_sympd()) {
        Rcout << "flag" << "\n";
        Ham_M.submat(2, 2, q+1, q+1) = X.t() * X;
    }
    
    arma::mat Ham_Minv = arma::inv_sympd(Ham_M, arma::inv_opts::allow_approx);
    int Ham_L = 4;
    double Ham_eps = 0.25;
    int MH_counter2 = 0;
    int TT_max = 5;
    
    arma::mat params_store(niter, q+2, arma::fill::zeros);
    arma::mat logpmat(10000, n, arma::fill::zeros);
    arma::vec sumloglik(10000, arma::fill::zeros);
    arma::vec ETA(n, arma::fill::zeros);
    arma::vec LAMBDA(n, arma::fill::zeros);
    arma::vec GAMMA(n, arma::fill::zeros);
    double ALPHA = 0.0;
    
    arma::vec PHI_init(q+2, arma::fill::zeros);
    arma::vec PHI(q+2, arma::fill::zeros);
    arma::vec params_ham(q+2, arma::fill::zeros);
    double logMH;
    double tempr;
    int tempiter;
    int i = 0;
    
    for (int iter = 0; iter < niter; iter++) {
        
        tempiter = iter + 1;
        
        // Estimate TT
        
        if ((tempiter - floor(tempiter/10)*10) == 0) {
            
            ETA = X * beta + reL + reR;
            for (unsigned int j = 0; j < ETA.n_rows; j++) {
                if (ETA[j] > 10) ETA[j] = 10;
            }
            LAMBDA = exp(ETA*(2-p))/(phi*(2-p));
            GAMMA = phi*(p-1)*exp((p-1)*ETA);
            ALPHA = (2-p)/(p-1);
            for (unsigned int i = 0; i < ids.size(); i++) {
                
                arma::vec TT_prob(TT_max, arma::fill::zeros);
                for (int j = 0; j < TT_max; j++) {
                    TT_prob[j] = gsl_ran_poisson_pdf((j+1), LAMBDA[ids[i]]) * gsl_ran_gamma_pdf(Y[ids[i]], ALPHA*(j+1), GAMMA[ids[i]]);
                }
                if(all(TT_prob == 0)) {
                    TT[ids[i]] = TT_max;
                } else {
                    arma::ivec temp = oneMultinomCalt(TT_prob);
                    TT[ids[i]] = temp.index_max() + 1;
                }
            }
            
        }
        
        // HMC step
        
        PHI_init = mvrnormArma_mean0(Ham_M);
        
        params_ham = params;
        arma::mat tempmat = Ham_eps * Ham_Minv;
        
        PHI = PHI_init + Ham_eps * logdens_dparams(params, X, Y, TT, beta_cov, reL, reR)*0.5;
        for (int L = 0; L < (Ham_L-1); L++) {
            params_ham = params_ham + tempmat * PHI;
            PHI = PHI + Ham_eps * logdens_dparams(params_ham, X, Y, TT, beta_cov, reL, reR);
        }
        params_ham = params_ham + tempmat * PHI;
        PHI = PHI + Ham_eps * logdens_dparams(params_ham, X, Y, TT, beta_cov, reL, reR)*0.5;
        
        if (params_ham.has_nan() || PHI.has_nan()) {
            logMH = datum::nan;
        } else {
            logMH = logdens(params_ham, X, Y, TT, beta_cov, reL, reR) - 
                logdens(params, X, Y, TT, beta_cov, reL, reR) +
                dmvnrm_diff(PHI, PHI_init, Ham_M, TRUE);
        }
        
        if (!isnan(logMH)) {
            tempr = gsl_ran_flat(s, 0, 1);
            if (tempr < exp(logMH)) {
                params = params_ham;
                MH_counter2 = MH_counter2 + 1;
                logphi = params[0];
                phi = exp(logphi);
                theta = params[1];
                if (theta > 10) {
                    theta = 10;
                } else if (theta < -10) {
                    theta = -10;
                }
                exptheta = exp(theta);
                p = (2*exptheta + 1) / (exptheta + 1);
                beta = params.rows(2, q+1);
            }
        }
        
        if( ((tempiter - floor(tempiter/100)*100) == 0) && (iter < 6501) ){
            
            // examine TT_max
            
            arma::uvec TT_check = find(TT.rows(ids) == TT_max);
            if( TT_check.size() > n*0.001 ) {
                TT_max = TT_max + 5;
            }
            
            // examine acceptance rate
            
            if( MH_counter2 > 85 ) {
                Ham_M = 0.65*0.65*Ham_M;
                Ham_Minv = arma::inv_sympd(Ham_M, arma::inv_opts::allow_approx);
            } else if( MH_counter2 > 75 ) {
                Ham_M = 0.75*0.75*Ham_M;
                Ham_Minv = arma::inv_sympd(Ham_M, arma::inv_opts::allow_approx);
            } else if( MH_counter2 > 65 ) {
                Ham_M = 0.85*0.85*Ham_M;
                Ham_Minv = arma::inv_sympd(Ham_M, arma::inv_opts::allow_approx);
            } else if( MH_counter2 < 25 ) {
                Ham_M = 1.3*1.3*Ham_M;
                Ham_Minv = arma::inv_sympd(Ham_M, arma::inv_opts::allow_approx);
            } else if( MH_counter2 < 35 ) {
                Ham_M = 1.2*1.2*Ham_M;
                Ham_Minv = arma::inv_sympd(Ham_M, arma::inv_opts::allow_approx);
            } else if( MH_counter2 < 45 ) {
                Ham_M = 1.1*1.1*Ham_M;
                Ham_Minv = arma::inv_sympd(Ham_M, arma::inv_opts::allow_approx);
            }
            MH_counter2 = 0;
        }
        
        if( (tempiter - floor(tempiter/500)*500) == 0 ){
            Rcout << tempiter << "\n";
            if( tempiter == 10000 ) {
                arma::vec params_diff = arma::diff(params_store.submat(1, 5, iter, 5));
                if( all(params_diff < 0.0001) ) break;
                if( MH_counter2 < 1600 ) {
                    break;
                } else {
                    MH_counter2 = 0; 
                }
            } else if( tempiter == 13000 ) {
                if(MH_counter2 < 1100) break;
            }
        }
        
        params_store.row(iter) = params.t();
        
        if(iter >= 10000) {
            
            arma::vec logpvec = loglik_ind(params, X, Y, TT, beta_cov, reL, reR);
            logpmat.row(i) = logpvec.t();
            sumloglik.row(i) = arma::accu(logpvec);
            i = i + 1;
            
        }
        
    }
    
    gsl_rng_free(s); // Free memory
    
    // for WAIC
    arma::rowvec avgloglik = arma::mean(logpmat, 0);
    arma::rowvec avglik = arma::mean(arma::exp(logpmat), 0);
    arma::rowvec varloglik = arma::var(logpmat, 0, 0);
    
    arma::vec waic(2, arma::fill::zeros);
    // WAIC 1
    waic[0] = 2 * arma::accu(log(avglik) - 2*avgloglik);
    // WAIC 2
    waic[1] = 2 * arma::accu(varloglik - log(avglik));
    
    return Rcpp::List::create(Rcpp::Named("params_store")=params_store,
                              Rcpp::Named("TT")=TT,
                              Rcpp::Named("lp")=sumloglik,
                              Rcpp::Named("WAIC")=waic);
    
}
