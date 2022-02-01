#ifndef __dwi_sdeconv_response_h__
#define __dwi_sdeconv_response_h__

#include "math/ZSH.h"
#include "math/math.h"
#include "mrtrix.h"

// Create a class to represent the response
//    public attributes of the class:
//      1. Method to retrieve the coeffs at a given b-value
//      2. Method to query the max harmonic order (that it stores)

// const MR::App::OptionGroup GradImportOptions = DWI::GradImportOptions();

namespace MR
{
  namespace DWI
  {
    namespace SDeconv
    {
      class Response
      {
        // add method to load matrixes and bvalues that correspond to it; + error handling ()

        // to make sure that bvals.size() = original.coefs.rows() to make sure the size is the same.
        //        - if it's not the case, we throw an exception: assert it and check for it

        // constructor or method taking a const str as arg; load the data from that;
        //        - create constructor with the same arguments as the method
        // get the comments, find the one with shells, strip the shells to get the nrs, parse_float for bvalues vector

      public:
        // template<typename VecType>
        Eigen::VectorXd retrieve_coeffs(const double bval)
        {
          Eigen::VectorXd coefs(original_coefs.cols());

          int i = 0;
          vector<double> dist_up, dist_down;

          while (bval < original_bvals[i])
            i++;

          // Linear interpolation
          float ratio = (bval - original_bvals[i]) / (original_bvals[i + 1] - original_bvals[i]);
          coefs = (1.0 - ratio) * original_coefs.row(i) + ratio * original_coefs.row(i + 1);

          if (i == original_coefs.cols())
            throw Exception("Iterator is out of boundaries");

          // if (!(original_coefs.coeffs().minCoeff() < bval < original_coefs.coeffs().maxCoeff()))
          //   throw Exception("b-value is out of boundaries, cannot interpolate");

          // // calculate difference between bval and elements of original_bvals
          // if (bval == i)
          //   return bval;
          // else if (i > bval)
          //   dist_up.push_back(abs(i - bval));
          // else
          //   dist_down.push_back(abs(i - bval));

          // // find bval neighbours
          // float up_bound = find(dist_up.begin(), dist_up.end(), min_element(dist_up.begin(), dist_up.end())),
          //       low_bound = find(dist_down.begin(), dist_down.end(), min_element(dist_down.begin(), dist_down.end()));

          // // interpolate
          // float ratio = (bval - low_bound) / (up_bound - low_bound);
          // coefs = (1.0 - ratio) * original_coefs.row(low_bound) + original_coefs.row(up_bound) * ratio;

          return coefs;
        }

        int lmax() const
        {

          return Math::ZSH::LforN(original_coefs.cols());
        }

        // LOAD MATRIX
        Eigen::MatrixXd load_matrix(const std::string &filename, vector<std::string>* comments)
        {
          auto vec = load_matrix_2D_vector(filename, comments);
          Eigen::MatrixXd M(vec.size(), vec[0].size());

          for (ssize_t r = 0; r < vec.size(); r++)
            for (ssize_t c = 0; c < vec[0].size(); r++)
              M(r, c) = vec[r][c];

          return M;
        }

        
        // const bool shell_bvalues = App::get_options("shell_bvalues").size();
        // const bool shell_sizes = App::get_options("shell_sizes").size();
        // const bool shell_indices = App::get_options("shell_indices").size();

        // Eigen::MatrixXd grad = DWI::get_DW_scheme(header, DWI::get_cmdline_bvalue_scaling_behaviour());

        // if (shell_bvalues || shell_sizes || shell_indices)
        //   print_shells(grad, shell_bvalues, shell_sizes, shell_indices);

      private:
        Eigen::MatrixXd original_coefs;
        vector<double> original_bvals;
      };
    }
  }
}

#endif
