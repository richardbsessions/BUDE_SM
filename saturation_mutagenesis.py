#!/usr/bin/env python3
# encoding: utf-8
'''
saturation_mutagenesis -- This utility takes a PDB file and Chain(s) ID for saturation mutagenesis

It expects BUDE(budeScan), BUDE sequence(bude_sequence) and Scwrl4 to be in the path.

It will mutate every residue in the chain(s) sequence to the list of residues 
to mutate to. Depending on the mode it will mutate the sequence's residues into
the full default list of mutable resiudes or those given in the command line. 

Sub commands

full:

Do a full mutagenesis with default mutable resiudes list.

manual:

Do mutagenesis with the given residues.

If the option to disable rotamer correction is given, the option to activate
residues for rotamer correction will be ignored.

@author:     Amaurys Ávila Ibarra

@copyright:  2019 University of Bristol. All rights reserved.

@license:    license

@contact:    Amaurys.AvilaIbarra@bristol.ac.uk
@deffield    updated: Updated
'''

import sys, os, argparse
from time import strftime, localtime

import ampal

from ampal_funcs.query_ampal import is_multi_model
from utils import config
from utils.check_executables import check_executables
from utils.do_mutagenesis import start_mutagenesis
from utils.f_system import init_directories

__all__ = []
__version__ = 0.1
__date__ = '2019-05-09'
__updated__ = '2019-05-09'
__sub_commands__ = "full|manual"


def main(argv=None):
    try:
        assert sys.version_info >= config.required_python
    except Exception as e:
        sys.stderr.write("\nFATAL ERROR: Wrong Python Version:\n\n"
            "We shall NOT continue.\nWe need python {}.{} or greater"
              " for using the budeAlaScan.\n".format(config.required_python[0],
                                                     config.required_python[1])
             )
        sys.stderr.write("\nWe found this Python Version:\n")
        sys.stderr.write("\nWe found this Python Version:\n")
        sys.stderr.write(sys.version)
        sys.stderr.write("\n")
        return 2

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_fulldesc = __import__('__main__').__doc__.split("\n")
    program_longdesc = '\n'.join(program_fulldesc[1:21])

    program_license = '''%s

  Created by Amaurys Ávila Ibarra on %s.
  Copyright 2019 University of Bristol. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.
   
   
 
USAGE
 
%s [%s] [options]

''' % (program_longdesc, str(__date__), program_name, __sub_commands__)

    full_desc = """
    This mode will do a full saturation mutagenesis.
    
    In this mode all residues in the sequence will be mutated to these default
    residues:
    [{}]
    
    If the option to disable rotamer correction is given, the option to activate
    residues for rotamer correction will be ignored.
""".format(", ".join(config.mutate_res))

    manual_desc = """
    This mode will do mutagenesis using the given residues.
    
    In this mode the mutagenesis will be done by mutating every residue in the
    sequence to the list of residues given in the command line.
    
    Resiudes given in the command line MUST be in this list.
    [{}]
    
    If the option to disable rotamer correction is given, the option to activate
    residues for rotamer correction will be ignored.
""".format(", ".join(config.legal_aa))
    try:
        # Setup argument parser
        parser = argparse.ArgumentParser(description=program_license, formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument('-v', '--verbose', dest="verbose_g", action='store_true', help="Output information during execution [default: %(default)s].")
        parser.add_argument('-t', '--turn-plot-off', dest="no_plot_g", action='store_false', help="Do not show plots [default: %(default)s].")
        parser.add_argument('-i', '--disable-rotamer-correction', dest="no_rot_g", action='store_true', help="No residues will be activated for rotamer correction [default: %(default)s].")

        subparsers = parser.add_subparsers(title="Mutagenesis sub-commands", description="Mutagenesis modes.", help="What mutagenesis mode to do.")

        # Full Mode
        parser_full = subparsers.add_parser('full', description=full_desc, help="Do mutagenesis with default list of mutable residues.", formatter_class=argparse.RawDescriptionHelpFormatter)
        parser_full.add_argument(dest="scan_mode", help=argparse.SUPPRESS, nargs='?', const=True, default='full')
        parser_full.add_argument('-v', '--verbose', dest="verbose", action='store_true', help="Output information during execution [default: %(default)s].")
        parser_full.add_argument('-t', '--turn-plot-off', dest="no_plot", action='store_false', help="Do not show plots [default: %(default)s].")

        # Mandatory options for full mode.
        mandatory_full = parser_full.add_argument_group("Mandatory Arguments")
        mandatory_full.add_argument("-p", "--pdb-file", dest="pdb_file_name", help="Name of the PDB file.", metavar="myPDB.pdb", required=True)
        mandatory_full.add_argument("-l", "--ligand-chains", dest="lig_chains", help="Ligand chain(s) to do the mutagenesis on.", metavar="A", type=str, nargs='+', required=True)

        # Default options for full mode.
        default_full = parser_full.add_argument_group("Default Arguments")
        default_full.add_argument("-a", "--residues-to-activate", dest="res2activate", help="Residues to activate for rotamer correction. [default: %(default)s].", metavar="D E", type=str, nargs='+', default="DERKH")
        default_full.add_argument('-i', '--disable-rotamer-correction', dest="no_rot_m", action='store_true', help="No residues will be activated for rotamer correction [default: %(default)s].")

        # manual selection
        parser_manual = subparsers.add_parser('manual', description=manual_desc, help="Do mutagenesis with the given resisues, [F M I L Y W].", formatter_class=argparse.RawDescriptionHelpFormatter)
        parser_manual.add_argument(dest="scan_mode", help=argparse.SUPPRESS, nargs='?', const=True, default='manual')
        parser_manual.add_argument('-v', '--verbose', dest="verbose", action='store_true', help="Output information during execution [default: %(default)s].")
        parser_manual.add_argument('-t', '--turn-plot-off', dest="no_plot", action='store_false', help="Do not show plots [default: %(default)s].")

        # Mandatory manual options.
        mandatory_manual = parser_manual.add_argument_group("Mandatory Arguments")
        mandatory_manual.add_argument("-p", "--pdb-file", dest="pdb_file_name", help="Name of the PDB file.", metavar="myPDB.pdb", required=True)
        mandatory_manual.add_argument("-l", "--ligand-chains", dest="lig_chains", help="Ligand chain(s) to do the mutagenesis on.", metavar="A", type=str, nargs='+', required=True)
        mandatory_manual.add_argument("-m", "--residues-for-mutagenesis", dest="mut_resisues", help="One letter code of residues to mutate to.", metavar="W", type=str, nargs='+', required=True)

        # Default options for full mode.
        default_manual = parser_manual.add_argument_group("Default Arguments")
        default_manual.add_argument("-a", "--residues-to-activate", dest="res2activate", help="Residues to activate for rotamer correction. [default: %(default)s].", metavar="D E", type=str, nargs='+', default="DERKH")
        default_manual.add_argument('-i', '--disable-rotamer-correction', dest="no_rot_m", action='store_true', help="No residues will be activated for rotamer correction [default: %(default)s].")
        # Process arguments
        args = parser.parse_args()

        config.pdb_basename = os.path.basename(args.pdb_file_name)[:-4]

        if args.verbose or args.verbose_g:
            config.verbose = True

        if args.no_rot_g or args.no_rot_m:
            config.do_rotamer_correction = False

        if not args.no_plot or not args.no_plot_g:
            config.showplots = False

        if config.verbose:
            print("{} started for PDB ID {} on: {}".format(program_name, config.pdb_basename, strftime(config.date_fmt, localtime())))
            if config.do_rotamer_correction:
                print("Rotamer correction is active.")
            else:
                print("Rotamer correction is NOT active.")

        if not check_executables():
            print("FATAL ERROR: Executable(s) not found.", file=sys.stderr)
            return 2
        if not os.path.isfile(args.pdb_file_name):
            print("FATAL ERROR: The PDB file [{}] cannot be found.".format(args.pdb_file_name), file=sys.stderr)
            return 2

        try:
            my_ampal = ampal.load_pdb(args.pdb_file_name)
        except Exception as e:
            sys.stderr.write("\nFatal Error: '{}' contains non standard amino acids.\n".format(args.pdb_file_name))
            sys.stderr.write("Error Info:\n{}\n".format(repr(e)))
            sys.exit(2)

        config.is_multimodel = is_multi_model(my_ampal)
        init_directories()

        start_mutagenesis(my_ampal, args)

        if config.verbose:
            print("{} Finished on: ".format(program_name), strftime(config.date_fmt, localtime()))

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help\n\n")
        return 2
    return 0


if __name__ == "__main__":
    # If no arguments is given show help.
    if len(sys.argv) is 1:
        sys.argv.append("-h")
    sys.exit(main())

