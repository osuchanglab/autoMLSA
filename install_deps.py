#!/usr/bin/env python
from util.validate_requirements import install_blast, install_mafft,\
    install_iqtree, init_logger, validate_requirements


def main():
    init_logger(False)
    install_blast()
    install_mafft()
    install_iqtree()
    validate_requirements()


if __name__ == '__main__':
    main()
