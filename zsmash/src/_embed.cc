#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <exception>
#include <boost/timer/timer.hpp>

#include <stdlib.h>
#include "config.h"
#include "semantic.h"


using namespace boost::program_options;

const string VERSION="sippl embed v1.0 \n Copyright I Chattopadhyay 2017";
const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";

vector<option> ignore_numbers(vector<string>& args)
{
  vector<option> result;
  int pos = 0;
  while(!args.empty())
    {
      const auto& arg = args[0];
      double num;
      if(boost::conversion::try_lexical_convert(arg, num))
	{
	  result.push_back(option());
	  option& opt = result.back();

	  opt.position_key = pos++;
	  opt.value.push_back(arg);
	  opt.original_tokens.push_back(arg);

	  args.erase(args.begin());
	}
      else
	break;
    }

  return result;
}

//####################################

int main(int argc, char *argv[])
{
  string hfilename="H.dat",  outfileE="outE.txt", outfileD="outD.txt";
  bool TIMER=true;

  options_description infor( "Program information");
  infor.add_options()
    ("help,h", "print help message.")
    ("version,V", "print version number");


  options_description usg( "Usage");
  usg.add_options()
    ("datafile,f",value<string>(), "data file")
    ("timer,t",value< bool >(), "display timer [1 (true)] ");

  options_description outputopt( "Output options");
  outputopt.add_options()
    ("Efile,E",value<string>(), "embedding coordinates file")
    ("Dfile,D",value<string>(), "dimension error file")
    ("verbose,v",value< int >(), "verbosity level [default: 1] ");

  options_description desc( "Sippl Embedding");
  desc.add(infor).add(usg).add(outputopt);

  positional_options_description p;
  variables_map vm;
  if (argc == 1)
    {
      cout << EMPTY_ARG_MESSAGE << endl;
      return 1;
    }
  try
    {
      store(command_line_parser(argc, argv)
	    .extra_style_parser(&ignore_numbers)
	    .options(desc)
	    .run(), vm);
      notify(vm);
    } 
  catch (std::exception &e)
    {
      cout << endl << e.what() 
	   << endl << desc << endl;
      return 1;
    }
  if (vm.count("help"))
    {
      cout << desc << endl;
      return 1;
    }
  if (vm.count("version"))
    {
      cout << VERSION << endl; 
      return 1;
    }

  if (vm.count("timer"))
    TIMER=vm["timer"].as<bool>();
  if (vm.count("datafile"))
    hfilename=vm["datafile"].as<string>();
  if (vm.count("Efile"))
    outfileE=vm["Efile"].as<string>();
  if (vm.count("Dfile"))
    outfileD=vm["Dfile"].as<string>();

  CONFIG cfg("config.cfg");


  matrix_dbl H;
  ifstream Hfile(hfilename.c_str());
  stringstream ss;

  string line;
  unsigned int row=0;
  while (getline (Hfile, line))
    {
      unsigned int col=0;
      map < unsigned int, double> Hr;
      double val;
      ss.clear();
      ss.str ("");
      ss << line;
      while ( ss>>val)
	Hr[col++]=val;
      H[row++]=Hr;
    }
  Hfile.close();
  
  if (TIMER)
    {
      timer::auto_cpu_timer t;

      Sippl_embed E(&H);
      matrix_dbl M2= E.embedding_coordinates(H.size());
      ofstream outE(outfileE.c_str());
      outE  << M2 << endl;
      outE.close();
      
      ofstream outD(outfileD.c_str());
      outD  << E.embedding_error() << endl;
      outD.close();
    }

  return 0;
}
