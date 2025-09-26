#include <regex>

#include <atrip/Checkpoint.hpp>

namespace atrip {

#define CHECKPOINT_FIELD_SEPARATION "\t";
namespace details {

/*
  ;;
  ;; regex created from the following recipee
  ;;
  (rxt-elisp-to-pcre
   (rx (or (and bol (0+ (any whitespace)) eol)
           (and bol (0+ any) (or "No" "Nv") (0+ any) eol)
           (and bol (0+ (any whitespace)) "It" (1+ any) eol))))
 */

const std::regex
    checkpoint_non_data_line("^[[:space:]]*$|^.*N[vo].*$|^[[:space:]]*It.+$");

} // namespace details

void write_checkpoint_header(std::string const &filepath) {
  std::ofstream out(filepath, std::ios::app | std::ios::out);
  for (auto const &name : std::vector<std::string>({"iteration",
                                                    "no",
                                                    "nv",
                                                    "nranks",
                                                    "nnodes",
                                                    "global_energy",
                                                    "iteration_energy",
                                                    "global_ct_energy",
                                                    "iteration_ct_energy",
                                                    "rank_round_robin"})) {
    out << name << CHECKPOINT_FIELD_SEPARATION;
  }
  // end newline
  out << "\n";
}

void write_checkpoint(Checkpoint const &c, std::string const &filepath) {
  const auto sep = CHECKPOINT_FIELD_SEPARATION;
  std::ofstream out(filepath, std::ios::app | std::ios::out);
  out                                                          /**/
      << c.iteration << sep                                    /**/
      << c.no << sep                                           /**/
      << c.nv << sep                                           /**/
      << c.nranks << sep                                       /**/
      << c.nnodes << sep                                       /**/
      << std::setprecision(19) << c.global_energy << sep       /**/
      << std::setprecision(19) << c.iteration_energy << sep    /**/
      << std::setprecision(19) << c.global_ct_energy << sep    /**/
      << std::setprecision(19) << c.iteration_ct_energy << sep /**/
      << c.rank_round_robin << sep                             /**/
      ;                                                        /**/
  // end newline
  out << "\n";
}

/**
 * @brief Reads checkpoint data from an input stream.
 *
 * This function reads checkpoint data line by line from the provided input
 * stream. It skips lines that match a predefined non-data line pattern. For
 * data lines, it extracts values for iteration, no, nv, nranks, nnodes,
 * global_energy, iteration_energy, and rank_round_robin, and populates a
 * Checkpoint object.
 *
 * @param in A reference to an input file stream from which to read the
 * checkpoint data.
 * @return A Checkpoint object containing the extracted data.
 */
Checkpoint read_checkpoint(std::ifstream &in) {
  Checkpoint c;
  for (std::string line; std::getline(in, line, '\n');) {
    // If it is a non data line then continue reading the file
    if (std::regex_match(line, details::checkpoint_non_data_line)) continue;
    std::stringstream s(line);
    s                            /**/
        >> c.iteration           /**/
        >> c.no                  /**/
        >> c.nv                  /**/
        >> c.nranks              /**/
        >> c.nnodes              /**/
        >> c.global_energy       /**/
        >> c.iteration_energy    /**/
        >> c.global_ct_energy    /**/
        >> c.iteration_ct_energy /**/
        >> c.rank_round_robin    /**/
        ;
  }
  return c;
}

Checkpoint read_checkpoint(std::string const &filepath) {
  std::ifstream in(filepath);
  return read_checkpoint(in);
}

} // namespace atrip
