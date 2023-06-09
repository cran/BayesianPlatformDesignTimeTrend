// Generated by rstantools.  Do not edit by hand.

/*
    BayesianPlatformDesignTimeTrend is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BayesianPlatformDesignTimeTrend is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BayesianPlatformDesignTimeTrend.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_logisticdummy_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_logisticdummy");
    reader.add_event(50, 48, "end", "model_logisticdummy");
    return reader;
}
#include <stan_meta_header.hpp>
class model_logisticdummy
  : public stan::model::model_base_crtp<model_logisticdummy> {
private:
        double beta0_prior_mu;
        double beta1_prior_mu;
        double beta0_prior_sigma;
        double beta1_prior_sigma;
        int beta0_nu;
        int beta1_nu;
        int N;
        std::vector<int> y;
        int K;
        std::vector<int> z;
        matrix_d x;
        int Kc;
        matrix_d xc;
        vector_d means_x;
public:
    model_logisticdummy(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_logisticdummy(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_logisticdummy_namespace::model_logisticdummy";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "beta0_prior_mu", "double", context__.to_vec());
            beta0_prior_mu = double(0);
            vals_r__ = context__.vals_r("beta0_prior_mu");
            pos__ = 0;
            beta0_prior_mu = vals_r__[pos__++];
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "beta1_prior_mu", "double", context__.to_vec());
            beta1_prior_mu = double(0);
            vals_r__ = context__.vals_r("beta1_prior_mu");
            pos__ = 0;
            beta1_prior_mu = vals_r__[pos__++];
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "beta0_prior_sigma", "double", context__.to_vec());
            beta0_prior_sigma = double(0);
            vals_r__ = context__.vals_r("beta0_prior_sigma");
            pos__ = 0;
            beta0_prior_sigma = vals_r__[pos__++];
            check_greater_or_equal(function__, "beta0_prior_sigma", beta0_prior_sigma, 0);
            current_statement_begin__ = 7;
            context__.validate_dims("data initialization", "beta1_prior_sigma", "double", context__.to_vec());
            beta1_prior_sigma = double(0);
            vals_r__ = context__.vals_r("beta1_prior_sigma");
            pos__ = 0;
            beta1_prior_sigma = vals_r__[pos__++];
            check_greater_or_equal(function__, "beta1_prior_sigma", beta1_prior_sigma, 0);
            current_statement_begin__ = 8;
            context__.validate_dims("data initialization", "beta0_nu", "int", context__.to_vec());
            beta0_nu = int(0);
            vals_i__ = context__.vals_i("beta0_nu");
            pos__ = 0;
            beta0_nu = vals_i__[pos__++];
            check_greater_or_equal(function__, "beta0_nu", beta0_nu, 0);
            current_statement_begin__ = 9;
            context__.validate_dims("data initialization", "beta1_nu", "int", context__.to_vec());
            beta1_nu = int(0);
            vals_i__ = context__.vals_i("beta1_nu");
            pos__ = 0;
            beta1_nu = vals_i__[pos__++];
            check_greater_or_equal(function__, "beta1_nu", beta1_nu, 0);
            current_statement_begin__ = 11;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 1);
            current_statement_begin__ = 12;
            validate_non_negative_index("y", "N", N);
            context__.validate_dims("data initialization", "y", "int", context__.to_vec(N));
            y = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("y");
            pos__ = 0;
            size_t y_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < y_k_0_max__; ++k_0__) {
                y[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 13;
            context__.validate_dims("data initialization", "K", "int", context__.to_vec());
            K = int(0);
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            K = vals_i__[pos__++];
            check_greater_or_equal(function__, "K", K, 1);
            current_statement_begin__ = 14;
            validate_non_negative_index("z", "N", N);
            context__.validate_dims("data initialization", "z", "int", context__.to_vec(N));
            z = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("z");
            pos__ = 0;
            size_t z_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < z_k_0_max__; ++k_0__) {
                z[k_0__] = vals_i__[pos__++];
            }
            size_t z_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < z_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "z[i_0__]", z[i_0__], 0);
                check_less_or_equal(function__, "z[i_0__]", z[i_0__], K);
            }
            current_statement_begin__ = 15;
            validate_non_negative_index("x", "N", N);
            validate_non_negative_index("x", "K", K);
            context__.validate_dims("data initialization", "x", "matrix_d", context__.to_vec(N,K));
            x = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, K);
            vals_r__ = context__.vals_r("x");
            pos__ = 0;
            size_t x_j_2_max__ = K;
            size_t x_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < x_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < x_j_1_max__; ++j_1__) {
                    x(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            // initialize transformed data variables
            current_statement_begin__ = 18;
            Kc = int(0);
            stan::math::fill(Kc, std::numeric_limits<int>::min());
            stan::math::assign(Kc,(K - 1));
            current_statement_begin__ = 19;
            validate_non_negative_index("xc", "N", N);
            validate_non_negative_index("xc", "Kc", Kc);
            xc = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, Kc);
            stan::math::fill(xc, DUMMY_VAR__);
            current_statement_begin__ = 20;
            validate_non_negative_index("means_x", "Kc", Kc);
            means_x = Eigen::Matrix<double, Eigen::Dynamic, 1>(Kc);
            stan::math::fill(means_x, DUMMY_VAR__);
            // execute transformed data statements
            current_statement_begin__ = 21;
            for (int i = 2; i <= K; ++i) {
                current_statement_begin__ = 22;
                stan::model::assign(means_x, 
                            stan::model::cons_list(stan::model::index_uni((i - 1)), stan::model::nil_index_list()), 
                            mean(stan::model::rvalue(x, stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list())), "x")), 
                            "assigning variable means_x");
                current_statement_begin__ = 23;
                stan::model::assign(xc, 
                            stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_uni((i - 1)), stan::model::nil_index_list())), 
                            subtract(stan::model::rvalue(x, stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list())), "x"), get_base1(means_x, (i - 1), "means_x", 1)), 
                            "assigning variable xc");
            }
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 27;
            validate_non_negative_index("beta1", "(K - 1)", (K - 1));
            num_params_r__ += (K - 1);
            current_statement_begin__ = 28;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_logisticdummy() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 27;
        if (!(context__.contains_r("beta1")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta1 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta1");
        pos__ = 0U;
        validate_non_negative_index("beta1", "(K - 1)", (K - 1));
        context__.validate_dims("parameter initialization", "beta1", "vector_d", context__.to_vec((K - 1)));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta1((K - 1));
        size_t beta1_j_1_max__ = (K - 1);
        for (size_t j_1__ = 0; j_1__ < beta1_j_1_max__; ++j_1__) {
            beta1(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta1);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta1: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 28;
        if (!(context__.contains_r("beta0")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta0 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta0");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "beta0", "double", context__.to_vec());
        double beta0(0);
        beta0 = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(beta0);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta0: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 27;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta1;
            (void) beta1;  // dummy to suppress unused var warning
            if (jacobian__)
                beta1 = in__.vector_constrain((K - 1), lp__);
            else
                beta1 = in__.vector_constrain((K - 1));
            current_statement_begin__ = 28;
            local_scalar_t__ beta0;
            (void) beta0;  // dummy to suppress unused var warning
            if (jacobian__)
                beta0 = in__.scalar_constrain(lp__);
            else
                beta0 = in__.scalar_constrain();
            // model body
            current_statement_begin__ = 35;
            lp_accum__.add(student_t_log(beta1, beta1_nu, beta1_prior_mu, beta1_prior_sigma));
            current_statement_begin__ = 36;
            lp_accum__.add(student_t_log(beta0, beta0_nu, beta0_prior_mu, beta0_prior_sigma));
            current_statement_begin__ = 39;
            lp_accum__.add(bernoulli_logit_glm_log(y, xc, beta0, beta1));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("beta1");
        names__.push_back("beta0");
        names__.push_back("b_Intercept");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back((K - 1));
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_logisticdummy_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta1 = in__.vector_constrain((K - 1));
        size_t beta1_j_1_max__ = (K - 1);
        for (size_t j_1__ = 0; j_1__ < beta1_j_1_max__; ++j_1__) {
            vars__.push_back(beta1(j_1__));
        }
        double beta0 = in__.scalar_constrain();
        vars__.push_back(beta0);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 44;
            double b_Intercept;
            (void) b_Intercept;  // dummy to suppress unused var warning
            stan::math::initialize(b_Intercept, DUMMY_VAR__);
            stan::math::fill(b_Intercept, DUMMY_VAR__);
            stan::math::assign(b_Intercept,(beta0 - dot_product(means_x, beta1)));
            // validate, write generated quantities
            current_statement_begin__ = 44;
            vars__.push_back(b_Intercept);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_logisticdummy";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta1_j_1_max__ = (K - 1);
        for (size_t j_1__ = 0; j_1__ < beta1_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta1" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "beta0";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "b_Intercept";
        param_names__.push_back(param_name_stream__.str());
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta1_j_1_max__ = (K - 1);
        for (size_t j_1__ = 0; j_1__ < beta1_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta1" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "beta0";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "b_Intercept";
        param_names__.push_back(param_name_stream__.str());
    }
}; // model
}  // namespace
typedef model_logisticdummy_namespace::model_logisticdummy stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
