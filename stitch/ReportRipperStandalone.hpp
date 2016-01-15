///////////////////////////////////////////////////////////////////////////////
// Copyright © 2013, Battelle Memorial Institute
//
// THIS FILE INITIALLY CREATED WITH:
//     TEMPLATE NAME:  lang_cpp_hdr.template
//     COMMAND NAME:   gensrc
//
// CODING CONVENTIONS:
//    * Class names are CamelCase with first word upper case
//    * Functions are camelCase with first word lower case
//    * Function parameters are lower case with _ and have p_ prefix
//    * Member variables always use 'this' pointer
///////////////////////////////////////////////////////////////////////////////


#ifndef REPORTRIPPERSTANDALONE_HPP
#define REPORTRIPPERSTANDALONE_HPP


// System Includes
#include <map>
#include <string>
#include <vector>
#include <cstdint>
// External Includes
// Internal Includes
// Application Includes


#ifdef  __CLASS__
#undef  __CLASS__
#endif
#define __CLASS__ "ReportRipperStandalone"


// Namespaces used
using namespace std;


typedef vector< pair<unsigned long int, unsigned long int> > hsp_type;
typedef struct hit_line
{
  unsigned long int ids;
  unsigned long int hsp_len;
  unsigned long int mismatch;
  unsigned long int gaps;
  unsigned long int q_start;
  unsigned long int q_end;
  unsigned long int s_start;
  unsigned long int s_end;
  unsigned long int merges;
  double    score;
  double    e_value;
  string    key;
  string    query;
  string    subject;
  hsp_type  q_hsps;
  hsp_type  s_hsps;
} hit_line;


typedef map<string, hit_line> HitLines;


class ReportRipperStandalone
{

public:

  // Getters & Setters

  inline HitLines getHitLines()
  {
    return(this->hit_lines);
  }


  // Constructors

  ReportRipperStandalone(string const & p_filename) throw();


  // Destructor

  ~ReportRipperStandalone(void) throw()
  {
    return;
  }


  // Public Methods

  string print(void) const throw();


private:

  // Getters/Setters


  // Constructors


  // Operators

  bool operator=(ReportRipperStandalone const & p_ripper) const throw();


  // Methods

  void merge(hit_line & p_existing, hit_line & p_tmp) throw();


  // Variables

  string   filename;
  HitLines hit_lines;


friend
  ostream & operator<<(ostream & p_os, ReportRipperStandalone const * p_ripper) throw();
friend
  ostream & operator<<(ostream & p_os, ReportRipperStandalone const & p_ripper) throw();

};




#endif // REPORTRIPPERSTANDALONE_HPP
