import pyopenmmpol as ommp

infile = "/home/mattia/LibEnv/open-mmpol/tests/N-methylacetamide/input_AMOEBA.mmp"
outfile = None

ommp.set_verbose(3)
ommp.init_mmp(infile)

if outfile is not None:
    ommp.print_summary_to_file(outfile)
else:
    ommp.print_summary()

ommp.terminate()
