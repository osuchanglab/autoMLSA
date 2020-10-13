#!/usr/bin/env python
from util.validate_requirements import install_blast, install_mafft,\
    install_iqtree, validate_requirements
from util.parse_args import init_logger


def main():
    init_logger(debug=False, quiet=False, rundir='', runid='')
    install_blast()
    install_mafft()
    install_iqtree()
    validate_requirements()


if __name__ == '__main__':
    main()
